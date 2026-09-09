"""
=======================================================================================
OT2 Protocol: ABTS Laccase Activity Assay - Tube-Rack Sourced (W-O-SC-02 / W-O-EN-01)
=======================================================================================
 Doc No:        W-O-SC-02 / W-O-EN-01
 Version:       v3.1
 Date:          2026-08-06
 Owner:         Lukas
 Workflow ID:   W-O-SC-02
 Unit ops:      U-O-01, U-O-02, U-O-03
 Output:        NUM_SAMPLES purified-laccase samples (in triplicate) + 1 shared
                no-enzyme control triplicate, each 200 uL, ready for the plate reader
                (A414 & A734). No mix/develop step runs on the robot -- the plate
                reader owns timing.

 SCOPE / DESIGN RATIONALE:
 - Two reagents only: purified laccase (10 uL/well) and ABTS reaction mix (190 uL/well,
   tops every well up to 200 uL).
 - Samples arrive pre-loaded in individual 1.5 mL tubes on the same style of
   3D-printed tube rack used by the GGA protocol (W-O-CL-01) -- pulled with the
   p20 in single-nozzle mode, one tube + one fresh tip per sample. Any dilution the
   samples need is prepared off-deck by the operator before the run.
 - One shared control triplicate (not per-sample): no-enzyme, i.e. water replaces the
   10 uL of enzyme, ABTS mix as normal. Doubles as the plate blank.
 - ABTS mix is added to whole columns via the p300 in full 8-channel mode --
   unused wells in a partially-filled last column also receive reagent; harmless
   since those wells are never sampled or read (same tradeoff already accepted in
   the GGA script for its last reaction column).
=======================================================================================
 STEPS TABLE
=======================================================================================
 U-O-01   Load enzyme from tubes into triplicate wells (P20, single-nozzle)
 U-O-02   Load water into the no-enzyme control wells (P20, single-nozzle)
 U-O-03   Add ABTS reaction mix to every well (P300, full 8-channel, reservoir)
  hand off to plate reader, operator reads A414 and A734 for x time
=======================================================================================
"""

import math
from opentrons import protocol_api
from opentrons.protocol_api import SINGLE

metadata = {
    'apiLevel': '2.20',
    'protocolName': 'Laccase - ABTS Activity Assay (Tube-Rack Sourced)',
    'description': (
        'Pulls purified laccase samples from individual tubes on a 3D-printed tube '
        'rack into triplicate assay-plate wells, adds one shared no-enzyme control '
        'triplicate, then tops every well up with ABTS reaction mix before handing '
        'off to the plate reader for A414 and A734.'
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

ENZYME_VOL_UL = 10.0
ABTS_MIX_VOL_UL = 190.0            # 10 + 190 = 200 uL final volume per well
WATER_VOL_UL = ENZYME_VOL_UL       # water replaces the enzyme in the no-enzyme control

# Plate layout: one triplicate per sample, plus the shared no-enzyme triplicate,
# filled column-then-row. The last column may be partially used; its unused wells
# still receive ABTS mix (see rationale above) but are never read.
NUM_GROUPS = NUM_SAMPLES + 1
NUM_WELLS_CONSUMED = NUM_GROUPS * 3
NUM_PLATE_COLUMNS = math.ceil(NUM_WELLS_CONSUMED / WELLS_PER_COLUMN)
PLATE_COLUMNS = [f'A{i + 1}' for i in range(NUM_PLATE_COLUMNS)]

# Tube rack (same style as the GGA script's, opentrons_24_tuberack_eppendorf_1.5ml_
# safelock_snapcap) has 4 ROWS (A-D), unlike the assay plate's 8 -- a separate row
# count avoids generating nonexistent well names like 'E1'.
TUBERACK_ROWS_PER_COLUMN = 4

DISPENSE_ASPIRATE_CLEARANCE_MM = 13.0

# Reservoir well assignments (nest_12_reservoir_15ml)
RESERVOIR_ABTS_WELL = 'A1'
RESERVOIR_WATER_WELL = 'A2'           # MQW for the no-enzyme control wells

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


# Consecutive triplicate groups: one per sample, then the no-enzyme control.
_names = _ordered_well_names(NUM_WELLS_CONSUMED)
_groups = [_names[i:i + 3] for i in range(0, len(_names), 3)]
SAMPLE_WELL_GROUPS = _groups[:NUM_SAMPLES]
NO_ENZYME_WELL_GROUP = _groups[NUM_SAMPLES]

TUBE_MAP = {_tube_well_name(j): SAMPLE_WELL_GROUPS[j] for j in range(NUM_SAMPLES)}


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
    # water sourcing only) -- configured once here, no later switch needed.
    # Anchored at 'H1' (front-most nozzle), not 'A1': the assay plate sits in slot 1
    # (a front-row slot), and 'A1' (the back-most nozzle) can't reach it without the
    # gantry exceeding its travel bounds -- PartialTipMovementNotAllowedError,
    # "outside of robot bounds" (same fix applied to the NNBT script's front-row
    # assay-plate access).
    p20.configure_nozzle_layout(style=SINGLE, start='H1', tip_racks=[tr_20])

    # Dispensed from DISPENSE_ASPIRATE_CLEARANCE_MM above the well bottom, not the
    # API's default ~1 mm clearance -- ABTS mix is added on top of the enzyme/water
    # already in each well, so dispensing close to that liquid risks the tip touching
    # it and cross-contaminating between columns.
    plate_columns = [
        assay_plate[w].bottom(z=DISPENSE_ASPIRATE_CLEARANCE_MM) for w in PLATE_COLUMNS]

    def _load_triplicate_single_nozzle(source_well, dest_well_names):
        """Pick up one tip, aspirate+dispense ENZYME_VOL_UL fresh from source_well for
        each of the 3 destination wells (the p20's 20 uL capacity can't hold all three
        doses at once), then drop the tip."""
        p20.pick_up_tip()
        for dest_name in dest_well_names:
            p20.aspirate(ENZYME_VOL_UL, source_well)
            p20.dispense(ENZYME_VOL_UL, assay_plate[dest_name])
        p20.drop_tip()

    # ===============================================================================
    # U-O-01  Load purified laccase from tubes (P20, single-nozzle, one tube + one
    #          reused tip per sample; fresh tip per sample -- different physical
    #          samples).
    #   Aspirates at the API's default clearance (no raised clearance) -- these are
    #   purified enzyme preps, no pelleted debris to stay clear of.
    # ===============================================================================
    for tube_well, dest_well_names in TUBE_MAP.items():
        _load_triplicate_single_nozzle(tuberack[tube_well], dest_well_names)

    # ===============================================================================
    # U-O-02  Load water into the no-enzyme control wells (P20, single-nozzle,
    #          reservoir source, one reused tip -- a shared triplicate, not per-sample)
    # ===============================================================================
    _load_triplicate_single_nozzle(reservoir[RESERVOIR_WATER_WELL], NO_ENZYME_WELL_GROUP)

    # ===============================================================================
    # U-O-03  Add ABTS reaction mix (P300, full 8-channel, whole columns) --
    #          reaction-start step.
    #   Tip reused across all columns -- same reagent to every destination, no
    #   cross-contamination concern.
    # ===============================================================================
    # One fresh aspirate per column, dispensed in full immediately -- no air_gap, no
    # distribute(). distribute() (and the air-gap/blow_out trick this replaced) loads
    # a single aspirate to serve several dispenses, and the air pocket between them
    # is exactly what let droplets cling to and slip off the tip between columns
    # (spillage/cross-contamination, same fix already applied to the NNBT script).
    # One aspirate matched to one dispense has no leftover air pocket to carry.
    # blow_out() is a separate call, not a dispense() keyword -- dispense() takes no
    # blow_out argument and raises TypeError if given one.
    p300.pick_up_tip()
    for column in plate_columns:
        p300.aspirate(ABTS_MIX_VOL_UL, reservoir[RESERVOIR_ABTS_WELL])
        p300.dispense(ABTS_MIX_VOL_UL, column)
        p300.blow_out(column)
    p300.drop_tip()
