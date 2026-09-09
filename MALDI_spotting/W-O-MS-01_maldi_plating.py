"""
=======================================================================================
OT2 Protocol: MALDI Target Plating - Matrix + Enzyme Product (W-O-MS-01)
=======================================================================================
 Doc No:        W-O-MS-01
 Version:       v1.1
 Date:          2026-09-08
 Owner:         Lukas
 Workflow ID:   W-O-MS-01
 Unit ops:      U-O-01, U-C-01
 Output:        `samples` enzyme-product samples spotted in `replicates` onto a
                maldi_384_wellplate target, each spot 1 uL matrix + 1 uL sample mixed
                in place, ready for the operator to dry and load into the MS.

 RUNTIME PARAMETERS (set in the Opentrons App before the run, no code edit needed):
   Samples            1-23,  default 3   number of enzyme-product tubes on the rack
   Start row          A-P,   default A   first target row to spot into
   Replicates         1-3,   default 1   spots per sample
   Matrix batch size  1-32,  default 9   samples per matrix->sample batch
 Everything else (volumes, mix strokes, z-height, deck slots) stays a constant in this
 file: those are bench-tuned values, not operator choices.

 SCOPE / DESIGN RATIONALE:
 - Samples and matrix both arrive pre-loaded in 1.5 mL tubes on the same style of
   3D-printed tube rack used by the GGA (W-O-CL-01) and ABTS (W-O-SC-02) protocols.
   The FIRST tube (A1) is the MALDI matrix; every tube after it is one enzyme-product
   sample. Any dilution the samples need is prepared off-deck before the run.
 - MATRIX DRIES. That is the whole reason this protocol is structured in batches
   rather than as two clean plate-wide passes (which is how the ABTS script is built).
   Spotting all replicates*samples matrix drops up front would leave the first ones
   sitting on the target for minutes before their sample lands. Instead the plate is
   worked in batches of `matrix_batch_size` samples: lay that batch's matrix spots,
   immediately add and mix that batch's samples, then move on. No matrix spot ever
   waits longer than one batch. The larger the batch, the longer the earliest spots in
   it sit -- the run log prints the drops-in-flight count so the operator can see it.
 - NO BLOW-OUT, NO AIR GAPS, ANYWHERE. Air pushed through a 1 uL droplet on a flat
   target sprays it across neighbouring spots and leaves bubbles in what is left, and
   a bubble in a 2 uL spot is a hole in the crystal layer. So: no blow_out(), no
   air_gap(), no distribute(), no touch_tip(). Every aspirate is matched 1:1 to a
   dispense, which is also what keeps a leftover air pocket from forming in the tip
   between drops (same reasoning as the ABTS script's rejection of distribute()).
 - The p20 is the only pipette on the deck. At 1 uL per drop it is working at its
   GEN2 minimum volume, which is also why the mix stroke below is 1.0 uL -- there is
   no legal smaller stroke, and a larger one would run the 2 uL spot dry mid-stroke
   and suck air back in.
=======================================================================================
 STEPS TABLE
=======================================================================================
 U-O-01   For each batch of `matrix_batch_size` samples:
            (a) spot 1 uL matrix into that batch's replicate wells
                (P20, single-nozzle, ONE tip for the whole batch)
            (b) for each sample in the batch, add 1 uL sample onto each of its matrix
                spots and mix 3x in place
                (P20, single-nozzle, FRESH tip per sample)
 U-C-01   PAUSE - hand off to the operator to dry the target and load the MS
=======================================================================================
"""

import math
from opentrons import protocol_api
from opentrons.protocol_api import SINGLE

metadata = {
    'protocolName': 'W-O-MS-01 - MALDI Target Plating (Matrix + Enzyme Product)',
    'description': (
        'Spots MALDI matrix and enzyme product from a 3D-printed tube rack onto a '
        '384-spot MALDI target, 1 uL + 1 uL per spot, mixed in place. Worked in small '
        'batches so no matrix spot dries before its sample arrives.'
    ),
}

# 2.18 is the minimum API level that supports runtime parameters; 2.20 is what the
# rest of the OT-2 protocols in this repo are on.
requirements = {'robotType': 'OT-2', 'apiLevel': '2.20'}

# =======================================================================================
# FIXED CONSTANTS  (bench-tuned -- deliberately NOT operator-settable)
# =======================================================================================
MATRIX_VOL_UL = 1.0
SAMPLE_VOL_UL = 1.0

# In-spot mixing, done by explicit aspirate/dispense rather than p20.mix() so the
# z-height of every stroke is under this file's control (mix() would otherwise work at
# whatever clearance the last motion left behind).
MIX_REPETITIONS = 3
MIX_VOL_UL = 1.0                   # the p20 GEN2's minimum legal volume, and half the
                                     # 2 uL sitting in the spot -- leaves 1 uL of
                                     # headroom so the stroke never bottoms out the
                                     # droplet and pulls air back in

# =======================================================================================
# Z-HEIGHT ON THE TARGET -- the number most likely to need bench tuning.
# The maldi_384_wellplate spots are 4 mm circles with a well DEPTH OF 0.1 mm: this is a
# flat target face, not a well. So a clearance measured from "well bottom" is very
# nearly a clearance measured from the bare plate surface, and the API's default ~1 mm
# would park the tip roughly 0.9 mm ABOVE the face -- clear of a 2 uL droplet, which
# would make the drop fall rather than be placed and would make the mix strokes
# aspirate air instead of liquid.
# 0.3 mm puts the tip ~0.2 mm off the face: into the droplet, but not pressed against
# the steel (which would occlude the tip orifice and dispense nothing).
# DRY-RUN THIS ON THE REAL TARGET BEFORE A REAL PLATING.
# =======================================================================================
MALDI_Z_CLEARANCE_MM = 0.3

# =======================================================================================
# TARGET LAYOUT -- maldi_384_wellplate, 16 rows (A-P) x 24 columns, 4.5 mm pitch.
# Replicate groups run LEFT TO RIGHT: sample 0 takes the first `replicates` columns of
# the start row, sample 1 the next `replicates`, and so on. 24 columns / replicates
# groups fill a row, after which the next sample wraps to column 1 of the row below.
# Rows above the start row are never touched.
# =======================================================================================
MALDI_ROW_LETTERS = [
    'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H',
    'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P']
MALDI_NUM_COLUMNS = 24

# Tube rack (same style as the GGA and ABTS scripts') has 4 ROWS (A-D), unlike the
# target's 16 -- a separate row count avoids generating nonexistent well names.
TUBERACK_ROWS_PER_COLUMN = 4
TUBERACK_NUM_COLUMNS = 6
TUBERACK_CAPACITY = TUBERACK_ROWS_PER_COLUMN * TUBERACK_NUM_COLUMNS   # 24 tubes
MATRIX_TUBE_INDEX = 0              # the FIRST tube on the rack holds the matrix;
                                     # samples occupy every tube after it

TIPS_PER_RACK = 96

# Deck slots. The same partial-tip-pickup geometry the ABTS script had to solve applies
# here: the target sits in slot 1 (front row) and is reached by the p20 in single-nozzle
# mode anchored at 'H1', which in turn keeps the p20's tip rack out of slot 4 and the
# tube rack off the slot adjacent to the tip rack. Verified with opentrons_simulate.
SLOT_MALDI_PLATE = '1'
SLOT_TUBERACK = '2'
SLOT_TIPRACK_20 = '7'

# The bench-agreed comfortable ceiling on matrix drops in flight at once. Above this the
# earliest matrix spots of a batch have been sitting for a while before their sample
# lands, so the protocol prints a warning into the run log rather than refusing -- the
# App's own limits (below) are what actually bound the parameter.
MATRIX_DROPS_IN_FLIGHT_ADVISORY = 9


# =======================================================================================
# RUNTIME PARAMETERS -- shown in the Opentrons App, set per run by the operator.
# =======================================================================================
def add_parameters(parameters: protocol_api.Parameters):

    parameters.add_int(
        variable_name='samples',
        display_name='Samples',
        description=(
            'Number of enzyme-product tubes on the rack (the matrix tube is extra).'),
        default=3,
        minimum=1,
        # The rack has 24 positions and the matrix tube takes one of them, so 23 is
        # the most that physically fits alongside it.
        maximum=TUBERACK_CAPACITY - 1,
    )

    parameters.add_str(
        variable_name='start_row',
        display_name='Start row',
        description=(
            'First target row to spot into, so a partly-used target can be continued.'),
        default='A',
        choices=[{'display_name': r, 'value': r} for r in MALDI_ROW_LETTERS],
    )

    parameters.add_int(
        variable_name='replicates',
        display_name='Replicates',
        description='Spots per sample, laid left to right in adjacent columns.',
        default=1,
        minimum=1,
        maximum=3,
    )

    parameters.add_int(
        variable_name='matrix_batch_size',
        display_name='Matrix batch size',
        description=(
            'Samples per batch. Bigger = fewer tips but longer matrix dry time.'),
        default=9,
        minimum=1,
        maximum=32,
    )


def _tube_well_name(index):
    """Tube-rack well name for the `index`-th physical tube (0-based) -- same
    row-then-column layout as the GGA script's GBLOCK_TUBE_MAP and the ABTS script's
    sample map (4 rows per column, A-D)."""
    col = index // TUBERACK_ROWS_PER_COLUMN
    row = MALDI_ROW_LETTERS[index % TUBERACK_ROWS_PER_COLUMN]
    return f'{row}{col + 1}'


def _spot_well_names(sample_index, start_row_index, replicates, groups_per_row):
    """The `replicates` target wells for the `sample_index`-th sample (0-based),
    running left to right from the start row and wrapping to the next row down every
    `groups_per_row` samples."""
    row_offset, position_in_row = divmod(sample_index, groups_per_row)
    row = MALDI_ROW_LETTERS[start_row_index + row_offset]
    first_column = position_in_row * replicates + 1
    return [f'{row}{first_column + i}' for i in range(replicates)]


# =======================================================================================
def run(protocol: protocol_api.ProtocolContext):

    num_samples = protocol.params.samples
    start_row = protocol.params.start_row
    replicates = protocol.params.replicates
    matrix_batch_size = protocol.params.matrix_batch_size

    # -----------------------------------------------------------------------------
    # DERIVED LAYOUT + VALIDATION
    #
    # These run during the App's protocol analysis, i.e. as soon as the operator picks
    # their parameter values and long before any tip moves -- so a combination that
    # would run off the bottom of the target or off the end of the rack is rejected on
    # screen, not mid-run on a half-spotted target. Ordered so each check runs BEFORE
    # the code that would otherwise blow up on the same bad value with a less useful
    # message.
    # -----------------------------------------------------------------------------
    start_row_index = MALDI_ROW_LETTERS.index(start_row)
    groups_per_row = MALDI_NUM_COLUMNS // replicates

    if MATRIX_TUBE_INDEX + 1 + num_samples > TUBERACK_CAPACITY:
        raise ValueError(
            f'{num_samples} samples + 1 matrix tube is more than the '
            f'{TUBERACK_CAPACITY} tubes on the rack. Split the run across two plates.')

    rows_needed = math.ceil(num_samples / groups_per_row)
    if start_row_index + rows_needed > len(MALDI_ROW_LETTERS):
        raise ValueError(
            f'{num_samples} samples at {replicates} replicates need {rows_needed} '
            f'rows, which runs off the bottom of the target starting from row '
            f'{start_row}. Start higher up, or reduce samples/replicates.')

    spot_volume_ul = MATRIX_VOL_UL + SAMPLE_VOL_UL
    if MIX_VOL_UL >= spot_volume_ul:
        raise ValueError(
            f'MIX_VOL_UL {MIX_VOL_UL} is the whole {spot_volume_ul} uL spot -- the '
            'stroke would bottom the droplet out and pull air back in. Leave headroom.')

    def spot_wells(sample_index):
        return _spot_well_names(
            sample_index, start_row_index, replicates, groups_per_row)

    sample_batches = [
        list(range(start, min(start + matrix_batch_size, num_samples)))
        for start in range(0, num_samples, matrix_batch_size)
    ]
    tips_needed = len(sample_batches) + num_samples
    if tips_needed > TIPS_PER_RACK:
        raise ValueError(
            f'{tips_needed} tips needed ({len(sample_batches)} matrix batches + '
            f'{num_samples} samples), more than the {TIPS_PER_RACK} on one rack.')

    # Sample i is drawn from the tube AFTER the matrix tube; the value is that sample's
    # group of target wells.
    tube_map = {
        _tube_well_name(MATRIX_TUBE_INDEX + 1 + i): spot_wells(i)
        for i in range(num_samples)
    }
    if len(tube_map) != num_samples:
        raise ValueError(
            'Sample tube assignments collided -- check samples against the rack size.')

    drops_in_flight = matrix_batch_size * replicates
    protocol.comment(
        f'Run plan: {num_samples} samples x {replicates} replicate(s) from row '
        f'{start_row}, in {len(sample_batches)} batch(es) of up to '
        f'{matrix_batch_size} samples ({drops_in_flight} matrix drops in flight), '
        f'{tips_needed} tips.')
    if drops_in_flight > MATRIX_DROPS_IN_FLIGHT_ADVISORY:
        protocol.comment(
            f'WARNING: {drops_in_flight} matrix drops per batch is above the '
            f'{MATRIX_DROPS_IN_FLIGHT_ADVISORY} checked on the bench -- the first '
            'spots of each batch may dry before their sample lands. Reduce the matrix '
            'batch size unless dry time has been re-checked.')

    # -----------------------------------------------------------------------------
    # DECK + INSTRUMENT SETUP
    # -----------------------------------------------------------------------------
    maldi_plate = protocol.load_labware('maldi_384_wellplate', SLOT_MALDI_PLATE)
    tuberack = protocol.load_labware('3d_printed_tuberack_1.5ml', SLOT_TUBERACK)  # TBD, see spec
    tr_20 = protocol.load_labware('opentrons_96_tiprack_20ul', SLOT_TIPRACK_20)

    p20 = protocol.load_instrument('p20_multi_gen2', 'right', tip_racks=[tr_20])

    # p20 stays in single-nozzle mode for the whole protocol. Anchored at 'H1'
    # (front-most nozzle), not 'A1': the target sits in slot 1 (a front-row slot), and
    # 'A1' (the back-most nozzle) cannot reach it without the gantry exceeding its
    # travel bounds -- PartialTipMovementNotAllowedError, "outside of robot bounds"
    # (same fix as the ABTS and NNBT scripts' front-row plate access).
    p20.configure_nozzle_layout(style=SINGLE, start='H1', tip_racks=[tr_20])

    matrix_tube = tuberack[_tube_well_name(MATRIX_TUBE_INDEX)]

    def _spot(well_name):
        """The point on the target every dispense and mix stroke works at -- see the
        MALDI_Z_CLEARANCE_MM comment for why this is not the API default."""
        return maldi_plate[well_name].bottom(z=MALDI_Z_CLEARANCE_MM)

    # ===============================================================================
    # U-O-01  Plate the target, one batch of `matrix_batch_size` samples at a time.
    #
    #   (a) MATRIX: one tip serves the whole batch -- it only ever sees the one matrix
    #       tube and clean target spots, so there is nothing to cross-contaminate. Each
    #       drop is a fresh 1 uL aspirate matched to one dispense; no multi-dispense
    #       from a single aspirate, because the air pocket that leaves behind is what
    #       makes droplets cling to and slip off the tip between spots.
    #
    #   (b) SAMPLE: fresh tip per sample -- these are different physical samples, and
    #       the tip enters the spot to mix, so it cannot be shared. Each of the
    #       sample's spots gets 1 uL dispensed straight into the matrix drop, then
    #       MIX_REPETITIONS strokes in place at the same low z.
    #
    #   Aspirates from the tubes use the API's default clearance -- these are clear
    #   liquids with no pellet to stay above.
    #
    #   Nothing here blows out, air-gaps or touches off. See the header.
    # ===============================================================================
    for batch_number, batch in enumerate(sample_batches, start=1):
        batch_wells = [w for i in batch for w in spot_wells(i)]

        protocol.comment(
            f'Batch {batch_number}/{len(sample_batches)}: matrix into '
            f"{', '.join(batch_wells)}")

        # (a) matrix for the whole batch, one tip
        p20.pick_up_tip()
        for well_name in batch_wells:
            p20.aspirate(MATRIX_VOL_UL, matrix_tube)
            p20.dispense(MATRIX_VOL_UL, _spot(well_name))
        p20.drop_tip()

        # (b) sample onto that matrix, fresh tip per sample, mixed in place
        for sample_index in batch:
            tube_well_name = _tube_well_name(MATRIX_TUBE_INDEX + 1 + sample_index)
            spot_names = tube_map[tube_well_name]

            protocol.comment(
                f'Batch {batch_number}/{len(sample_batches)}: sample from tube '
                f"{tube_well_name} into {', '.join(spot_names)}")

            p20.pick_up_tip()
            for well_name in spot_names:
                target = _spot(well_name)
                p20.aspirate(SAMPLE_VOL_UL, tuberack[tube_well_name])
                p20.dispense(SAMPLE_VOL_UL, target)
                for _ in range(MIX_REPETITIONS):
                    p20.aspirate(MIX_VOL_UL, target)
                    p20.dispense(MIX_VOL_UL, target)
            p20.drop_tip()

    # ===============================================================================
    # U-C-01  PAUSE - hand off to the operator.
    # ===============================================================================
    protocol.pause(
        f'MALDI target plated: {num_samples} samples x {replicates} replicate(s), '
        f'rows {start_row} onward. Let the spots dry, then load the target into the MS.')
