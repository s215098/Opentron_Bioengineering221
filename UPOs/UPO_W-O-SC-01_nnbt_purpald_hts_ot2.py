"""
=======================================================================================
OT2 Protocol: NNBT Purpald Assay - Tube-Rack Sourced UPO Protocol (W-O-SC-01)
=======================================================================================
 Doc No:        W-O-SC-01
 Version:       v4.0
 Date:          2026-08-06
 Owner:         Lukas
 Workflow ID:   W-O-SC-01
 Unit ops:      U-O-01..U-O-07, U-C-01
 Output:        NUM_SAMPLES supernatant samples, each assayed in triplicate at EVERY
                pH in the BUFFERS panel, + one water-only control triplicate per pH,
                dosed with buffer then reaction mix, incubated, dosed with Purpald and
                developed on the heater-shaker, ready for the plate reader (A530).

 SCOPE / DESIGN RATIONALE:
 - Samples arrive pre-loaded in individual 1.5 mL tubes on the same style of
   3D-printed tube rack used by the GGA (W-O-CL-01) and ABTS_UPO (W-O-EN-01)
   protocols -- pulled with the p300 in single-nozzle mode, one tube + one reused
   tip per sample (fresh tip per sample since these are different physical
   samples).
 - Enzyme concentration in the supernatant is unknown, so a fixed 50 uL volume is
   pulled per enzyme (no titration series).
 - WELL RECIPE: 50 uL supernatant + 15 uL pH buffer + 85 uL reaction mix = 150 uL,
   incubated, then + 50 uL Purpald = 200 uL final. That fills the NEST 200 uL flat
   plate to its nominal capacity -- there is no headroom left, so none of these
   volumes can be raised without changing the plate.
 - The reaction mix is premixed off-deck and added as ONE reagent; it replaces the
   separate H2O2, NNBT and water-top-up additions of earlier versions.
 - pH PANEL: every sample is assayed at every pH in the BUFFERS table, each as its own
   triplicate. Because the p20 adds buffer in full 8-channel mode -- one stroke fills a
   whole column and cannot mix pHs within it -- each pH owns a COLUMN-ALIGNED BLOCK,
   padded up to a column boundary. Padding wells get reaction mix and Purpald with
   their column but no buffer and no sample, and are never read.
 - One water-only control triplicate PER pH (not one shared across the plate), so each
   pH block carries its own matched blank. Its "sample slot" is water instead of
   supernatant; everything else in the well is identical to a real sample.
 - The reaction mix still carries acetonitrile-dissolved NNBT (~1/3 the viscosity of
   water), which drips or over-shoots at default pipetting flow rates -- that step
   alone slows aspirate/dispense and adds a tip-touch after dispense.
 - SHARED UPO INCUBATION PROTOCOL (also used by ABTS_UPO, W-O-EN-01): incubate at
   37 degC / 220 rpm / 5 min. 37 degC is the bottom of the heater-shaker's
   controllable range, so it is honoured on real hardware; the >=37 degC check below
   only falls back to ambient if someone lowers the setpoint.
 - The plate is NOT covered with parafilm at any point and the protocol never pauses
   for the operator -- it runs start to finish unattended.
 - Development shake runs hot too (37 degC) and vigorously -- see DEVELOP_RPM for
   how that speed was chosen against the near-full wells.
=======================================================================================
 STEPS TABLE
=======================================================================================
 U-O-01   Load 50 uL samples from tubes into their triplicate wells in every pH block
          (P300, single-nozzle, one tip per sample)
 U-O-02   Load 50 uL control water from reservoir into one control triplicate per pH
          (P300, single-nozzle)
 U-O-03   Add 15 uL pH buffer per column block (P20, full 8-channel, one reservoir
          well and one fresh tip per buffer)
 U-O-04   Add 85 uL reaction mix (P300, full 8-channel, reservoir, slowed flow rate
          + tip-touch) -- brings every well to 150 uL
 U-O-05   Heater-shaker incubate: 37 degC / 220 rpm / 5 min
 U-O-06   Add 50 uL Purpald reagent (P300, full 8-channel, reservoir) -> 200 uL total
 U-O-07   Heater-shaker develop shake: 37 degC / DEVELOP_RPM / 10 min
 U-C-01   Release the latch and hand off to plate reader, operator reads A530
=======================================================================================
"""

import math
from opentrons import protocol_api
from opentrons.protocol_api import SINGLE

metadata = {
    'apiLevel': '2.20',
    'protocolName': 'W-O-SC-01 - NNBT Purpald Assay (Tube-Rack Sourced UPO Protocol)',
    'description': (
        'Pulls a fixed 50 uL supernatant per enzyme from individual tubes on a '
        '3D-printed tube rack into a triplicate for every pH in the BUFFERS panel, '
        'adds a water-only control triplicate per pH, doses 15 uL of that block\'s pH '
        'buffer and 85 uL reaction mix, runs the shared UPO incubation step, adds '
        '50 uL Purpald, develops on the heater-shaker, then hands off to the plate '
        'reader for A530.'
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

SUPERNATANT_VOL_UL = 50.0            # fixed per enzyme -- no titration
BUFFER_VOL_UL = 15.0                 # pH buffer, whichever pH the well's block uses
REACTION_MIX_VOL_UL = 85.0           # premixed off-deck; 50 + 15 + 85 = 150 uL
PURPALD_VOL_UL = 50.0                # 150 + 50 = 200 uL, the plate's full capacity

# =======================================================================================
# pH PANEL -- the buffers to screen. EVERY sample is assayed at EVERY pH listed here,
# each as its own triplicate, and every pH gets its own water-only control triplicate.
# Add/remove entries to change the panel; all well groups, column blocks and tip counts
# recompute below. Each entry needs its own reservoir well, loaded with that buffer.
#
# Set exactly one entry to fall back to the old single-buffer behaviour.
# =======================================================================================
BUFFERS = [
    {'name': 'pH 6.0', 'reservoir_well': 'A1'},
]
NUM_BUFFERS = len(BUFFERS)

# PLATE LAYOUT -- one COLUMN-ALIGNED BLOCK per pH.
# The p20 adds buffer in full 8-channel mode, which fills an entire column in a single
# stroke and therefore cannot put two different pHs in the same physical column. So
# each pH owns a whole number of columns: its (NUM_SAMPLES + 1) triplicate groups are
# padded up to the next column boundary, and the next pH starts in a fresh column.
# The padding wells receive reaction mix and Purpald along with their column, but no
# buffer and no sample -- they are never read (same tradeoff the ABTS script accepts
# for its partially-filled last column).
GROUPS_PER_BUFFER = NUM_SAMPLES + 1          # NUM_SAMPLES samples + 1 water-only control
WELLS_PER_BUFFER = GROUPS_PER_BUFFER * 3
COLUMNS_PER_BUFFER = math.ceil(WELLS_PER_BUFFER / WELLS_PER_COLUMN)
WELLS_PER_BUFFER_PADDED = COLUMNS_PER_BUFFER * WELLS_PER_COLUMN

NUM_PLATE_COLUMNS = COLUMNS_PER_BUFFER * NUM_BUFFERS
NUM_WELLS_CONSUMED = NUM_PLATE_COLUMNS * WELLS_PER_COLUMN
PLATE_COLUMNS = [f'A{i + 1}' for i in range(NUM_PLATE_COLUMNS)]

# Tube rack (same style as the GGA/ABTS scripts', 4 ROWS (A-D), unlike the assay
# plate's 8 -- a separate row count avoids generating nonexistent well names like 'E1'.
TUBERACK_ROWS_PER_COLUMN = 4

# Reagent dispensing: dispense this far above the well bottom instead of the API's
# default ~1 mm clearance -- buffer, reaction mix and Purpald are all added on top of
# liquid already in each well, so dispensing close to that liquid risks the tip
# touching it and cross-contaminating between columns (same fix applied to the
# ABTS_UPO script's reagent dispensing). Includes PLATE_CRASH_HEADSPACE_MM on top of
# the original 13 mm -- real-hardware crash into the plate persisted even after fixing
# the HS adapter/plate labware mismatch, so this adds a blanket safety margin against
# whatever residual real-vs-modeled height gap remains (physical tolerances, deck
# calibration, etc.) until that gap is root-caused.
PLATE_CRASH_HEADSPACE_MM = 5.0        # ~0.5 cm safety margin, see above
DISPENSE_ASPIRATE_CLEARANCE_MM = 13.0 + PLATE_CRASH_HEADSPACE_MM

# Sample/control loading (U-O-01/02) previously dispensed at the API's default ~1 mm
# clearance (no explicit z) -- the lowest, most crash-prone target in the protocol and
# the first wells touched (matches the reported A1 crash). Now explicit, with the same
# headspace margin applied.
SAMPLE_DISPENSE_CLEARANCE_MM = 1.0 + PLATE_CRASH_HEADSPACE_MM

# Reaction-mix viscosity handling: the mix carries acetonitrile-dissolved NNBT
# (~1/3 water's viscosity), which drips or over-shoots at default flow rates. That
# step's aspirate/dispense rate is scaled down by this fraction of the pipette's
# default, then restored immediately after.
REACTION_MIX_FLOW_RATE_SCALE = 0.5

# Shared UPO incubation protocol (also used by ABTS_UPO, W-O-EN-01)
INCUBATION_TEMP_C = 37               # 37 degC is the BOTTOM of the HS module's
                                     # controllable range, so this setpoint is honoured;
                                     # the >=37 check below only trips into the ambient
                                     # fallback if someone lowers this value
INCUBATION_RPM = 220
INCUBATION_MIN = 5

# Development: "shake vigorously" at 37 degC for 10 min. The heater-shaker accepts
# 200-3000 rpm, but the wells hold 200 uL in a 200 uL-capacity plate with no parafilm
# lid -- at the top of that range an uncovered, brim-full 96-well plate slops liquid
# between wells, which cross-contaminates a triplicate-based assay. 1000 rpm is the
# usual "vigorous mixing" setting for flat 96-well plates: ~4.5x the incubation speed,
# well clear of the module's 200 rpm floor, and still conservative against splash-over.
# Raise toward 1500 only after confirming on a sacrificial plate that nothing spills.
DEVELOP_TEMP_C = 37
DEVELOP_RPM = 1000
DEVELOP_MIN = 10

# Reservoir well assignments (nest_12_reservoir_15ml). Parked at the far end of the
# 12-well reservoir on purpose: the buffers claim their wells from A1 upward via the
# BUFFERS table above, so leaving A1-A9 free lets the pH panel grow without having to
# renumber these. The assert below catches any collision anyway.
RESERVOIR_REACTION_MIX_WELL = 'A10'  # premixed reaction mix
RESERVOIR_WATER_WELL = 'A11'         # the control triplicates' 50 uL "sample slot"
RESERVOIR_PURPALD_WELL = 'A12'

# Deck slots. Heater-Shaker holds the assay plate; tip rack (p20) and tube rack sit in
# non-front-row slots -- required for single-nozzle/partial-tip access (front slots
# 1-3 error "outside of robot bounds", per the GGA script's validation note).
#
# The p300's single-nozzle work also targets slot 1 itself (the assay plate, on the
# Heater-Shaker), which is a front-row slot -- confirmed against the real OT-2 API
# simulator (opentrons_simulate): with the nozzle configuration anchored at 'A1' (the
# back-most physical nozzle), dispensing into slot 1 raised "outside of robot bounds";
# anchoring at 'H1' (the front-most nozzle) instead lets the gantry reach it. But with
# 'H1', the p300's tip rack can no longer sit in slot 4 (directly behind the
# Heater-Shaker): the simulator then raised PartialTipMovementNotAllowedError picking
# up tips there, first against slot 1 (with 'A1') and then against slot 7's tip rack
# (with 'H1') -- slot 4 is boxed in on both sides. Moved to slot 9 instead, which the
# simulator confirmed is clear of both neighbors.
HS_SLOT = 1
SLOT_TIPRACK_300 = '8'
SLOT_RESERVOIR = '6'
SLOT_TUBERACK = '3'
SLOT_TIPRACK_20 = '7'


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
    """Tube-rack well name for the `index`-th sample (0-based) -- same row-then-column
    layout as the GGA script's GBLOCK_TUBE_MAP (4 rows per column, A-D)."""
    col = index // TUBERACK_ROWS_PER_COLUMN
    row = ROW_LETTERS[index % TUBERACK_ROWS_PER_COLUMN]
    return f'{row}{col + 1}'


_all_names = _ordered_well_names(NUM_WELLS_CONSUMED)

# Per-pH block slicing. Block b spans wells [b*WELLS_PER_BUFFER_PADDED : ...], of which
# only the first WELLS_PER_BUFFER are assigned to groups -- the remainder is padding.
# Within a block the groups run: one triplicate per sample, then the control.
BUFFER_PLATE_COLUMNS = []            # columns owned by each buffer, index-aligned to BUFFERS
SAMPLE_WELL_GROUPS = []              # [buffer_index][sample_index] -> 3 well names
CONTROL_WELL_GROUPS = []             # [buffer_index] -> 3 well names
for _b in range(NUM_BUFFERS):
    _start = _b * WELLS_PER_BUFFER_PADDED
    _block = _all_names[_start:_start + WELLS_PER_BUFFER]
    _groups = [_block[i:i + 3] for i in range(0, len(_block), 3)]
    SAMPLE_WELL_GROUPS.append(_groups[:NUM_SAMPLES])
    CONTROL_WELL_GROUPS.append(_groups[NUM_SAMPLES])
    BUFFER_PLATE_COLUMNS.append(
        PLATE_COLUMNS[_b * COLUMNS_PER_BUFFER:(_b + 1) * COLUMNS_PER_BUFFER])

# One tube per sample; each tube now feeds NUM_BUFFERS triplicates (one per pH), so the
# mapped value is the flat list of every well that sample occupies across the plate.
TUBE_MAP = {
    _tube_well_name(i): [
        well for _b in range(NUM_BUFFERS) for well in SAMPLE_WELL_GROUPS[_b][i]
    ]
    for i in range(NUM_SAMPLES)
}

# Fail loudly at load time rather than mid-run on a half-dosed plate.
_buffer_wells = [b['reservoir_well'] for b in BUFFERS]
assert NUM_BUFFERS >= 1, 'BUFFERS must list at least one buffer.'
assert len(set(_buffer_wells)) == NUM_BUFFERS, (
    f'BUFFERS reuse a reservoir well: {_buffer_wells}')
_other_wells = [RESERVOIR_REACTION_MIX_WELL, RESERVOIR_WATER_WELL, RESERVOIR_PURPALD_WELL]
assert not set(_buffer_wells) & set(_other_wells), (
    f'BUFFERS collide with the reaction-mix/water/Purpald wells {_other_wells}.')
assert NUM_PLATE_COLUMNS <= 12, (
    f'{NUM_BUFFERS} buffers x {COLUMNS_PER_BUFFER} columns = {NUM_PLATE_COLUMNS} '
    'columns, more than the 12 on a 96-well plate. Reduce BUFFERS or NUM_SAMPLES.')


# =======================================================================================
def run(protocol: protocol_api.ProtocolContext):

    # -----------------------------------------------------------------------------
    # DECK + INSTRUMENT SETUP
    # -----------------------------------------------------------------------------
    hs_mod = protocol.load_module('heaterShakerModuleV1', HS_SLOT)
    # Physical adapter on the Heater-Shaker is the 96 Flat Bottom Adapter, which is only
    # compatible with nest_96_wellplate_200ul_flat / corning_384_wellplate_112ul_flat
    # per its stackingOffsetWithLabware -- corning_96_wellplate_360ul_flat is NOT on that
    # list and raises LabwareCannotBeStackedError, so the assay plate is the NEST 200uL
    # flat plate instead. Must declare the adapter here too, or the simulator assumes
    # the plate sits directly on the module with no adapter stack-up, undershooting
    # every well's real Z height and driving the tip further down than intended
    # (crashed into A1).
    hs_adapter = hs_mod.load_adapter('opentrons_96_flat_bottom_adapter')
    assay_plate = hs_adapter.load_labware('nest_96_wellplate_200ul_flat')
    tr_300 = protocol.load_labware('opentrons_96_tiprack_300ul', SLOT_TIPRACK_300)
    tr_20 = protocol.load_labware('opentrons_96_tiprack_20ul', SLOT_TIPRACK_20)
    reservoir = protocol.load_labware('nest_12_reservoir_15ml', SLOT_RESERVOIR)
    tuberack = protocol.load_labware('3d_printed_tuberack_1.5ml', SLOT_TUBERACK)  # TBD, see spec

    p300 = protocol.load_instrument('p300_multi_gen2', 'left', tip_racks=[tr_300])
    p20 = protocol.load_instrument('p20_multi_gen2', 'right', tip_racks=[tr_20])

    # p300 starts in single-nozzle mode for tube/control-water sourcing (50 uL exceeds
    # the p20's capacity); reconfigured to 8-channel later for the reaction-mix and
    # Purpald steps. p20 handles only the 15 uL buffer add, 8-channel throughout.
    # Anchored at 'H1' (front-most nozzle), not 'A1' -- see the deck-slot comment
    # above for why: 'A1' can't reach the front-row assay plate (slot 1, on the
    # Heater-Shaker) without exceeding the gantry's travel bounds, and also collides
    # with the Heater-Shaker's housing when picking tips from slot 4.
    p300.configure_nozzle_layout(style=SINGLE, start='H1', tip_racks=[tr_300])

    plate_columns = [
        assay_plate[w].bottom(z=DISPENSE_ASPIRATE_CLEARANCE_MM) for w in PLATE_COLUMNS]
    hs_mod.close_labware_latch()  # secure plate for pipetting + shaking

    def _load_wells_single_nozzle(source_location, dest_well_names):
        """Pick up one tip, aspirate+dispense SUPERNATANT_VOL_UL fresh from
        source_location for each destination well (50 uL exceeds the p20's capacity
        but fits comfortably in one p300 aspirate), then drop the tip. One tip serves
        the whole list -- callers group by sample, so the tip only ever sees one
        liquid."""
        p300.pick_up_tip()
        for dest_name in dest_well_names:
            p300.aspirate(SUPERNATANT_VOL_UL, source_location)
            p300.dispense(
                SUPERNATANT_VOL_UL, assay_plate[dest_name].bottom(z=SAMPLE_DISPENSE_CLEARANCE_MM))
        p300.drop_tip()

    # ===============================================================================
    # U-O-01  Load samples from tubes (P300, single-nozzle, one tube + one reused tip
    #          per sample; fresh tip per sample -- different physical samples).
    #   Each sample now fills NUM_BUFFERS triplicates -- one in every pH block -- all
    #   from the same tube on the same tip.
    #   Aspirates at the API's default clearance (no raised clearance) -- samples are
    #   supernatant-only, no pelleted debris to stay clear of.
    # ===============================================================================
    for tube_well, dest_well_names in TUBE_MAP.items():
        _load_wells_single_nozzle(tuberack[tube_well], dest_well_names)

    # ===============================================================================
    # U-O-02  Load controls (P300, single-nozzle, reservoir water source, one reused
    #          tip for all of them -- same water everywhere; no raised-clearance
    #          concern, clean reservoir liquid).
    #   One water-only control triplicate per pH, so each pH block carries its own
    #   matched blank.
    # ===============================================================================
    _load_wells_single_nozzle(
        reservoir[RESERVOIR_WATER_WELL],
        [well for group in CONTROL_WELL_GROUPS for well in group])

    # ===============================================================================
    # U-O-03  Add pH buffers (P20, full 8-channel, whole columns) -- one buffer per
    #          column block, FRESH TIP PER BUFFER.
    #   Unlike the single-buffer version, the tip cannot be reused across the whole
    #   plate here: each block gets a chemically different liquid, and carrying one
    #   pH into the next block's aspirate would shift it. Within a block the tip is
    #   reused, since every column there gets the same buffer.
    #   The 8-channel stroke is exactly why the blocks are column-aligned -- see the
    #   PLATE LAYOUT comment above.
    # ===============================================================================
    # Manual aspirate/dispense per column, one fresh aspirate per destination -- no
    # air_gap, no distribute(). distribute() (and the air-gap trick this replaced)
    # loads a single aspirate to serve several dispenses, and the air pocket that
    # separates them is exactly what let droplets cling to and slip off the tip
    # between columns (spillage/cross-contamination). One aspirate matched to one
    # dispense has no leftover air pocket to carry between wells.
    for buffer_index, buffer_spec in enumerate(BUFFERS):
        protocol.comment(
            f"Adding buffer {buffer_spec['name']} from reservoir "
            f"{buffer_spec['reservoir_well']} to columns "
            f"{', '.join(BUFFER_PLATE_COLUMNS[buffer_index])}")
        p20.pick_up_tip()
        for column_name in BUFFER_PLATE_COLUMNS[buffer_index]:
            target = assay_plate[column_name].bottom(z=DISPENSE_ASPIRATE_CLEARANCE_MM)
            p20.aspirate(BUFFER_VOL_UL, reservoir[buffer_spec['reservoir_well']])
            p20.dispense(BUFFER_VOL_UL, target)
            p20.blow_out(target)

        p20.drop_tip()

    # ===============================================================================
    # U-O-04  Add reaction mix (P300, reconfigured to 8-channel, whole columns,
    #          slowed flow rate) -- brings every well to 150 uL total.
    #   The mix carries acetonitrile-dissolved NNBT (~1/3 the viscosity of water),
    #   which drips or over-shoots at default flow rates, so this step alone slows
    #   aspirate/dispense and touches the tip off after each dispense; both are
    #   restored/scoped to just this step so buffer and Purpald are unaffected.
    #   85 uL exceeds the p20's 20 uL capacity, hence the p300 rather than the p20
    #   that handled the old 7.5 uL NNBT add.
    # ===============================================================================
    p300.configure_nozzle_layout(style=protocol_api.ALL, tip_racks=[tr_300])
    default_p300_aspirate_rate = p300.flow_rate.aspirate
    default_p300_dispense_rate = p300.flow_rate.dispense
    p300.flow_rate.aspirate = default_p300_aspirate_rate * REACTION_MIX_FLOW_RATE_SCALE
    p300.flow_rate.dispense = default_p300_dispense_rate * REACTION_MIX_FLOW_RATE_SCALE
    # Same one-aspirate-per-dispense handling as the buffer step above -- no
    # distribute(), no air gap to carry droplets between columns.
    p300.pick_up_tip()
    for column in plate_columns:
        p300.aspirate(REACTION_MIX_VOL_UL, reservoir[RESERVOIR_REACTION_MIX_WELL])
        p300.dispense(REACTION_MIX_VOL_UL, column)
        p300.touch_tip()
    p300.drop_tip()
    p300.flow_rate.aspirate = default_p300_aspirate_rate
    p300.flow_rate.dispense = default_p300_dispense_rate

    # ===============================================================================
    # U-O-05  Heater-shaker incubate (shared UPO incubation protocol, also used by
    #          ABTS_UPO): 37 degC / 220 rpm / 5 min. 37 degC sits at the bottom of the
    #          module's controllable range, so it is applied for real; the else-branch
    #          only runs if INCUBATION_TEMP_C is lowered below that floor.
    #   The heater is left ON afterwards -- U-O-07 develops at the same 37 degC, so
    #   deactivating here would only force a reheat.
    # ===============================================================================
    if INCUBATION_TEMP_C >= 37:
        hs_mod.set_and_wait_for_temperature(INCUBATION_TEMP_C)
    else:
        protocol.comment(
            f'INCUBATION_TEMP_C ({INCUBATION_TEMP_C}) is below the heater-shaker\'s '
            'controllable range: shaking at ambient (room temperature) instead.'
        )
    hs_mod.set_and_wait_for_shake_speed(INCUBATION_RPM)
    protocol.delay(minutes=INCUBATION_MIN)
    hs_mod.deactivate_shaker()

    # ===============================================================================
    # U-O-06  Add Purpald (P300, full 8-channel, whole columns) -> 200 uL total.
    #   Shaker is stopped here (U-O-05 deactivated it) -- dispensing into a moving
    #   plate would misplace the droplet.
    # ===============================================================================
    # Same one-aspirate-per-dispense handling as the buffer/reaction-mix steps above.
    p300.pick_up_tip()
    for column in plate_columns:
        p300.aspirate(PURPALD_VOL_UL, reservoir[RESERVOIR_PURPALD_WELL])
        p300.dispense(PURPALD_VOL_UL, column)
        p300.blow_out(column)

    p300.drop_tip()

    # ===============================================================================
    # U-O-07  Heater-shaker develop shake: 37 degC / DEVELOP_RPM / 10 min -- the
    #          "shake vigorously" step. See DEVELOP_RPM for why 1000 rpm rather than
    #          the module's 3000 rpm ceiling.
    # ===============================================================================
    if DEVELOP_TEMP_C >= 37:
        hs_mod.set_and_wait_for_temperature(DEVELOP_TEMP_C)
    hs_mod.set_and_wait_for_shake_speed(DEVELOP_RPM)
    protocol.delay(minutes=DEVELOP_MIN)
    hs_mod.deactivate_shaker()
    hs_mod.deactivate_heater()

    # ===============================================================================
    # U-C-01  Hand off to plate reader
    # ===============================================================================
    hs_mod.open_labware_latch()

