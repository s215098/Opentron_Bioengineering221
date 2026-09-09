"""
=======================================================================================
OT2 Protocol: ABTS_UPO Activity Assay - Tube-Rack Sourced (W-O-SC-02 / W-O-EN-01)
=======================================================================================
 Doc No:        W-O-SC-02 / W-O-EN-01
 Version:       v3.0
 Date:          2026-08-07
 Owner:         Lukas
 Workflow ID:   W-O-SC-02
 Unit ops:      U-O-01, U-O-02, U-O-03, U-O-04, U-C-01
 Output:        NUM_SAMPLES supernatant samples, each assayed in triplicate against
                EVERY ABTS mix in the ABTS_MIXES panel, + one no-enzyme control
                triplicate per mix, dosed with that block's ABTS mix and H2O2, ready
                for the plate reader (A414). No mix/develop step runs on the robot --
                the plate reader owns timing.

 SCOPE / DESIGN RATIONALE:
 - Samples arrive pre-loaded in individual 1.5 mL tubes on the same style of
   3D-printed tube rack used by the GGA protocol (W-O-CL-01) -- pulled with the
   p20 in single-nozzle mode, one tube + one reused tip per sample (fresh tip per
   sample since these are different physical samples). Any dilution the samples need
   is prepared off-deck by the operator before the run.
 - BUFFER PANEL: unlike the NNBT assay (W-O-SC-01), the buffer here is not a separate
   addition -- it is already inside the ABTS reaction mix. So a pH panel is run by
   loading several ABTS mixes, each premade in a different buffer, into their own
   reservoir wells (see ABTS_MIXES below). Every sample is assayed against every mix,
   each as its own triplicate. The 20 + 160 + 20 uL well recipe is untouched by this.
 - Because the p300 adds ABTS mix in full 8-channel mode -- one stroke fills a whole
   column and cannot mix two buffers within it -- each mix owns a COLUMN-ALIGNED
   BLOCK, padded up to a column boundary. Padding wells receive ABTS mix and H2O2
   with their column but no sample, and are never read (same tradeoff already
   accepted in the GGA and NNBT scripts).
 - One no-enzyme control triplicate PER mix (not one shared across the plate), so each
   buffer carries its own matched blank: water replaces the supernatant, everything
   else in the well is identical to a real sample. This is the same single control the
   NNBT script uses; the earlier no-H2O2 and blank triplicates are gone.
 - Every well on the plate therefore receives real H2O2, which is what lets U-O-04 be
   one blanket 8-channel pass. The earlier version had to split the plate into
   H2O2-recipient and water-instead column blocks precisely because the no-H2O2/blank
   controls existed; with those dropped, that split is unnecessary.
=======================================================================================
 STEPS TABLE
=======================================================================================
 U-O-01   Load samples from tubes into their triplicate wells in every mix block
          (P20, single-nozzle, one tip per sample)
 U-O-02   Load water into one no-enzyme control triplicate per mix (P20,
          single-nozzle, reservoir)
 U-O-03   Add ABTS reaction mix per column block (P300, full 8-channel, one reservoir
          well and one fresh tip per mix)
 U-O-04   Add H2O2 to every column (P300, full 8-channel, reservoir)
 U-C-01   PAUSE - hand off to plate reader, operator reads A414
=======================================================================================
"""

import math
from opentrons import protocol_api
from opentrons.protocol_api import SINGLE

metadata = {
    'apiLevel': '2.20',
    'protocolName': 'ABTS_UPO - ABTS Activity Assay (Tube-Rack Sourced)',
    'description': (
        'Pulls supernatant samples from individual tubes on a 3D-printed tube rack '
        'into a triplicate for every ABTS mix in the ABTS_MIXES panel, adds a '
        'no-enzyme control triplicate per mix, doses that block\'s ABTS mix and H2O2, '
        'then hands off to the plate reader for A414.'
    ),
}

# =======================================================================================
# PARAMETERS  (single source of truth - no magic numbers anywhere else in this file)
# =======================================================================================
NUM_SAMPLES = 3                    # placeholder -- single source of truth for batch
                                     # size; change this and TUBE_MAP / well groups /
                                     # column counts all recompute automatically below
WELLS_PER_COLUMN = 8
ROW_LETTERS = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']

SAMPLE_VOL_UL = 20.0
ABTS_MIX_VOL_UL = 160.0
H2O2_VOL_UL = 20.0
WATER_VOL_UL = SAMPLE_VOL_UL          # same 20 uL, standing in for the supernatant in
                                       # the no-enzyme control wells

# =======================================================================================
# BUFFER PANEL -- the ABTS mixes to screen. The buffer is part of the reaction mix in
# this assay (there is no separate buffer addition, unlike the NNBT script), so a pH
# panel means several ABTS mixes, each premade off-deck in a different buffer and loaded
# into its own reservoir well. EVERY sample is assayed against EVERY mix listed here,
# each as its own triplicate, and every mix gets its own no-enzyme control triplicate.
# Add/remove entries to change the panel; all well groups, column blocks and tip counts
# recompute below.
#
# Set exactly one entry to fall back to the old single-mix behaviour.
# =======================================================================================
ABTS_MIXES = [
    {'name': 'pH 4.0', 'reservoir_well': 'A1'},
]
NUM_MIXES = len(ABTS_MIXES)

# PLATE LAYOUT -- one COLUMN-ALIGNED BLOCK per ABTS mix.
# The p300 adds ABTS mix in full 8-channel mode, which fills an entire column in a
# single stroke and therefore cannot put two different buffers in the same physical
# column. So each mix owns a whole number of columns: its (NUM_SAMPLES + 1) triplicate
# groups are padded up to the next column boundary, and the next mix starts in a fresh
# column. The padding wells receive ABTS mix and H2O2 along with their column, but no
# sample -- they are never read.
GROUPS_PER_MIX = NUM_SAMPLES + 1             # NUM_SAMPLES samples + 1 no-enzyme control
WELLS_PER_MIX = GROUPS_PER_MIX * 3
COLUMNS_PER_MIX = math.ceil(WELLS_PER_MIX / WELLS_PER_COLUMN)
WELLS_PER_MIX_PADDED = COLUMNS_PER_MIX * WELLS_PER_COLUMN

NUM_PLATE_COLUMNS = COLUMNS_PER_MIX * NUM_MIXES
NUM_WELLS_CONSUMED = NUM_PLATE_COLUMNS * WELLS_PER_COLUMN
PLATE_COLUMNS = [f'A{i + 1}' for i in range(NUM_PLATE_COLUMNS)]

# Single-nozzle dispensing chunk: the p20's 20 uL capacity can't hold all three 20 uL
# doses for one sample/control at once, so each tip re-aspirates before every dispense.
SINGLE_NOZZLE_CHUNK_UL = SAMPLE_VOL_UL

# Tube rack (same style as the GGA script's, opentrons_24_tuberack_eppendorf_1.5ml_
# safelock_snapcap) has 4 ROWS (A-D), unlike the assay plate's 8 -- a separate row
# count avoids generating nonexistent well names like 'E1'.
TUBERACK_ROWS_PER_COLUMN = 4

DISPENSE_ASPIRATE_CLEARANCE_MM = 13.0


# Reservoir well assignments (nest_12_reservoir_15ml). Parked at the far end of the
# 12-well reservoir on purpose: the ABTS mixes claim their wells from A1 upward via the
# ABTS_MIXES table above, so leaving A1-A10 free lets the panel grow without having to
# renumber these. The asserts below catch any collision anyway.
RESERVOIR_H2O2_WELL = 'A11'
RESERVOIR_WATER_WELL = 'A12'          # the no-enzyme control triplicates' "sample slot"

# Deck slots. All verified against the real OT-2 API simulator (opentrons_simulate),
# not just reasoned about -- partial-tip-pickup deck restrictions turned out to be
# stricter, and different, than plain aspirate/dispense-from-a-well restrictions:
#
# - The assay plate (slot 1) is a front-row slot the p20 must reach in single-nozzle
#   mode. Anchored at 'A1' (back-most nozzle), dispensing there raised "outside of
#   robot bounds"; anchored at 'H1' (front-most nozzle) instead, it works.
# - With 'H1', the p20's own TIP RACK could no longer sit in slot 4 (directly next to
#   slot 1) -- PartialTipMovementNotAllowedError there. Moved to slot 9, confirmed
#   clear.
# - But slot 9 then boxed in the TUBE RACK in slot 6 (directly adjacent) the same way
#   -- PartialTipMovementNotAllowedError picking up tube contents there. Moved the
#   tube rack to slot 2 (front row) instead: unlike tip pickup, plain aspirate/
#   dispense from tube-rack wells has no front-row restriction, confirmed by a full
#   clean simulation run with the tube rack there.
SLOT_ASSAY_PLATE = '1'
SLOT_TIPRACK_20 = '7'
SLOT_RESERVOIR = '6'
SLOT_TUBERACK = '2'
SLOT_TIPRACK_300 = '8'


def _ordered_well_names(num_wells):
    """First `num_wells` well names of the assay plate, column-then-row order
    (A1, B1, ... H1, A2, B2, ...)."""
    names = []
    for col in range(NUM_PLATE_COLUMNS):
        for row in ROW_LETTERS:
            names.append(f'{row}{col + 1}')
            if len(names) == num_wells:
                return names
    return names


def _tube_well_name(index):
    """Tube-rack well name for the `index`-th physical tube (0-based) -- same
    row-then-column layout as the GGA script's GBLOCK_TUBE_MAP (4 rows per column,
    A-D)."""
    col = index // TUBERACK_ROWS_PER_COLUMN
    row = ROW_LETTERS[index % TUBERACK_ROWS_PER_COLUMN]
    return f'{row}{col + 1}'


_all_names = _ordered_well_names(NUM_WELLS_CONSUMED)

# Per-mix block slicing. Block m spans wells [m*WELLS_PER_MIX_PADDED : ...], of which
# only the first WELLS_PER_MIX are assigned to groups -- the remainder is padding.
# Within a block the groups run: one triplicate per sample, then the no-enzyme control.
MIX_PLATE_COLUMNS = []               # columns owned by each mix, index-aligned to ABTS_MIXES
SAMPLE_WELL_GROUPS = []              # [mix_index][sample_index] -> 3 well names
CONTROL_WELL_GROUPS = []             # [mix_index] -> 3 well names
for _m in range(NUM_MIXES):
    _start = _m * WELLS_PER_MIX_PADDED
    _block = _all_names[_start:_start + WELLS_PER_MIX]
    _groups = [_block[i:i + 3] for i in range(0, len(_block), 3)]
    SAMPLE_WELL_GROUPS.append(_groups[:NUM_SAMPLES])
    CONTROL_WELL_GROUPS.append(_groups[NUM_SAMPLES])
    MIX_PLATE_COLUMNS.append(
        PLATE_COLUMNS[_m * COLUMNS_PER_MIX:(_m + 1) * COLUMNS_PER_MIX])

# One tube per sample; each tube now feeds NUM_MIXES triplicates (one per buffer), so
# the mapped value is the flat list of every well that sample occupies across the plate.
TUBE_MAP = {
    _tube_well_name(i): [
        well for _m in range(NUM_MIXES) for well in SAMPLE_WELL_GROUPS[_m][i]
    ]
    for i in range(NUM_SAMPLES)
}

# Fail loudly at load time rather than mid-run on a half-dosed plate.
_mix_wells = [m['reservoir_well'] for m in ABTS_MIXES]
assert NUM_MIXES >= 1, 'ABTS_MIXES must list at least one mix.'
assert len(set(_mix_wells)) == NUM_MIXES, (
    f'ABTS_MIXES reuse a reservoir well: {_mix_wells}')
_other_wells = [RESERVOIR_H2O2_WELL, RESERVOIR_WATER_WELL]
assert not set(_mix_wells) & set(_other_wells), (
    f'ABTS_MIXES collide with the H2O2/water wells {_other_wells}.')
assert NUM_PLATE_COLUMNS <= 12, (
    f'{NUM_MIXES} mixes x {COLUMNS_PER_MIX} columns = {NUM_PLATE_COLUMNS} columns, '
    'more than the 12 on a 96-well plate. Reduce ABTS_MIXES or NUM_SAMPLES.')


# =======================================================================================
def run(protocol: protocol_api.ProtocolContext):

    # -----------------------------------------------------------------------------
    # DECK + INSTRUMENT SETUP
    # -----------------------------------------------------------------------------
    assay_plate = protocol.load_labware('corning_96_wellplate_360ul_flat', SLOT_ASSAY_PLATE)
    tr_20 = protocol.load_labware('opentrons_96_tiprack_20ul', SLOT_TIPRACK_20)
    tr_300 = protocol.load_labware('opentrons_96_tiprack_300ul', SLOT_TIPRACK_300)
    reservoir = protocol.load_labware('nest_12_reservoir_15ml', SLOT_RESERVOIR)
    tuberack = protocol.load_labware('3d_printed_tuberack_1.5ml', SLOT_TUBERACK)  # TBD, see spec

    p20 = protocol.load_instrument('p20_multi_gen2', 'right', tip_racks=[tr_20])
    p300 = protocol.load_instrument('p300_multi_gen2', 'left', tip_racks=[tr_300])

    # p20 stays in single-nozzle mode for its entire role in this protocol (tube +
    # control-water sourcing only) -- configured once here, no later switch needed.
    # Anchored at 'H1' (front-most nozzle), not 'A1': the assay plate sits in slot 1
    # (a front-row slot), and 'A1' (the back-most nozzle) can't reach it without the
    # gantry exceeding its travel bounds -- PartialTipMovementNotAllowedError,
    # "outside of robot bounds" (same fix applied to the NNBT script's front-row
    # assay-plate access).
    p20.configure_nozzle_layout(style=SINGLE, start='H1', tip_racks=[tr_20])

    # Dispensed from DISPENSE_ASPIRATE_CLEARANCE_MM above the well bottom, not the
    # API's default ~1 mm clearance -- ABTS mix and H2O2 are added on top of the
    # sample/water already in each well, so dispensing close to that liquid risks the
    # tip touching it and cross-contaminating between columns.
    plate_columns = [
        assay_plate[w].bottom(z=DISPENSE_ASPIRATE_CLEARANCE_MM) for w in PLATE_COLUMNS]

    def _load_wells_single_nozzle(source_well, dest_well_names):
        """Pick up one tip, aspirate+dispense SAMPLE_VOL_UL fresh from source_well for
        each destination well (the p20's 20 uL capacity can't hold two doses at once),
        then drop the tip. One tip serves the whole list -- callers group by sample, so
        the tip only ever sees one liquid."""
        p20.pick_up_tip()
        for dest_name in dest_well_names:
            p20.aspirate(SINGLE_NOZZLE_CHUNK_UL, source_well)
            p20.dispense(SINGLE_NOZZLE_CHUNK_UL, assay_plate[dest_name])
        p20.drop_tip()

    # ===============================================================================
    # U-O-01  Load samples from tubes (P20, single-nozzle, one tube + one reused tip
    #          per sample; fresh tip per sample -- different physical samples).
    #   Each sample now fills NUM_MIXES triplicates -- one in every buffer block -- all
    #   from the same tube on the same tip.
    #   Aspirates at the API's default clearance (no raised clearance) -- samples are
    #   supernatant-only, no pelleted debris to stay clear of.
    # ===============================================================================
    for tube_well, dest_well_names in TUBE_MAP.items():
        _load_wells_single_nozzle(tuberack[tube_well], dest_well_names)

    # ===============================================================================
    # U-O-02  Load the no-enzyme controls (P20, single-nozzle, reservoir water source,
    #          one reused tip for all of them -- same water everywhere).
    #   One control triplicate per mix, so each buffer block carries its own matched
    #   blank.
    # ===============================================================================
    _load_wells_single_nozzle(
        reservoir[RESERVOIR_WATER_WELL],
        [well for group in CONTROL_WELL_GROUPS for well in group])

    # ===============================================================================
    # U-O-03  Add ABTS reaction mix (P300, full 8-channel, whole columns) -- one mix
    #          per column block, FRESH TIP PER MIX.
    #   Unlike the single-mix version, the tip cannot be reused across the whole plate
    #   here: each block gets a chemically different liquid, and carrying one buffer
    #   into the next block's aspirate would shift its pH. Within a block the tip is
    #   reused, since every column there gets the same mix.
    #   The 8-channel stroke is exactly why the blocks are column-aligned -- see the
    #   PLATE LAYOUT comment above.
    # ===============================================================================
    # One fresh aspirate per column, dispensed in full immediately -- no air_gap, no
    # distribute(). distribute() (and the air-gap/blow_out trick this replaced) loads
    # a single aspirate to serve several dispenses, and the air pocket between them
    # is exactly what let droplets cling to and slip off the tip between columns
    # (spillage/cross-contamination, same fix already applied to the NNBT script).
    # One aspirate matched to one dispense has no leftover air pocket to carry.
    # blow_out() is a separate call, not a dispense() keyword -- dispense() takes no
    # blow_out argument and raises TypeError if given one.
    for mix_index, mix_spec in enumerate(ABTS_MIXES):
        protocol.comment(
            f"Adding ABTS mix {mix_spec['name']} from reservoir "
            f"{mix_spec['reservoir_well']} to columns "
            f"{', '.join(MIX_PLATE_COLUMNS[mix_index])}")
        p300.pick_up_tip()
        for column_name in MIX_PLATE_COLUMNS[mix_index]:
            target = assay_plate[column_name].bottom(z=DISPENSE_ASPIRATE_CLEARANCE_MM)
            p300.aspirate(ABTS_MIX_VOL_UL, reservoir[mix_spec['reservoir_well']])
            p300.dispense(ABTS_MIX_VOL_UL, target)
            p300.blow_out(target)
        p300.drop_tip()

    # ===============================================================================
    # U-O-04  Add H2O2 to every column (P300, full 8-channel) -- reaction-start step.
    #   Blanket, not column-scoped: the only control is the no-enzyme one, which gets
    #   real H2O2 like every sample well, so no column needs a different reagent here.
    #   (The earlier no-H2O2/blank controls are what forced the old split into
    #   H2O2-recipient and water-instead column blocks.)
    #   Tip reused across all columns -- same reagent to every destination, no
    #   cross-contamination concern.
    # ===============================================================================
    # Same one-aspirate-per-dispense handling as the ABTS mix step above.
    p300.pick_up_tip()
    for column in plate_columns:
        p300.aspirate(H2O2_VOL_UL, reservoir[RESERVOIR_H2O2_WELL])
        p300.dispense(H2O2_VOL_UL, column)
        p300.blow_out(column)
    p300.drop_tip()
