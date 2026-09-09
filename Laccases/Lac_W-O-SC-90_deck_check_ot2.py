"""
=======================================================================================
OT2 Protocol: DECK CHECK - hardware dry run for W-O-SC-03 / W-O-SC-04
=======================================================================================
 Doc No:        W-O-SC-90
 Version:       v1.0
 Owner:         Kristine Toft Johansen s215098

 WHAT THIS IS
 ------------
 A short, deliberately dumb protocol that makes the robot perform EVERY KIND OF MOVE
 the real assay needs - once each, with water - so that a geometry or calibration
 problem shows up in five minutes instead of halfway through a two-hour run with real
 enzyme on the deck.

 It is NOT an assay. It measures nothing. Load water everywhere.

 It is intentionally NOT built on W-O-SC-03's layout engine: it is short enough to read
 top to bottom, and its job is to test the hardware, not the plate layout.

 WHAT IT CHECKS, in increasing order of risk
 -------------------------------------------
   1  labware + module load, adapter stack-up, latch                 (no movement)
   2  P20 single-nozzle -> assay plate A1 on the Heater-Shaker       (the historic crash)
   3  P20 single-nozzle <- tube rack
   4  P20 single-nozzle <- dilution plate
   5  P20 fill + MIX in a dilution well   (regression test: mixing must
      draw liquid, not air - the tip has under 2 mm of depth to work in)
   6  P300 8-channel -> assay plate whole column  (p300 is only ever 8-channel)
   7  Heater-Shaker heat + both shake speeds
   8  latch release

 HOW TO USE IT
 -------------
   1. Fill reservoir A3 with ~5 mL water. Put ~200 uL water in tube rack A1.
      Everything else can be empty and dry.
   2. Put the SAME labware in the SAME slots you will use tomorrow. That is the point.
   3. opentrons_simulate -L Labware Lac_W-O-SC-90_deck_check_ot2.py
   4. Run it on the robot and WATCH. It pauses before each phase so you can put your
      hand on the E-stop and abort if anything looks low.
   5. Anything that fails here would have failed tomorrow. Fix it now.

 If a phase crashes or looks too low, set that phase's flag to False to skip past it
 and keep testing the rest, then deal with the failure separately.
=======================================================================================
"""

import math

from opentrons import protocol_api
from opentrons.protocol_api import SINGLE, ALL

metadata = {
    'apiLevel': '2.20',
    'protocolName': 'W-O-SC-90 - deck check (water only, no assay)',
    'description': (
        'Dry run for the NNBT/Purpald workflow: exercises every pipetting move and the '
        'Heater-Shaker once each, with water, pausing before each phase so the operator '
        'can watch. Not an assay.'
    ),
}

# ---------------------------------------------------------------------------------
# Phases. Set any of these False to skip that phase - useful for getting past a known
# failure to test whatever comes after it.
# ---------------------------------------------------------------------------------
CHECK_P20_TO_PLATE = True        # the historic A1 crash
CHECK_P20_FROM_TUBERACK = True
CHECK_P20_FROM_DILUTION = True
CHECK_P300_TO_DILUTION = True    # p20 multi-stroke fill + mix
CHECK_P300_MULTI_COLUMN = True
CHECK_HEATER_SHAKER = True

# Pause and wait for the operator before every phase. Leave True the first time.
STOP_BETWEEN_PHASES = True

# ---------------------------------------------------------------------------------
# These MUST match W-O-SC-03 exactly, or this check proves nothing about tomorrow.
# If you change a slot or a clearance there, change it here too.
# ---------------------------------------------------------------------------------
HS_SLOT = 1
SLOT_TUBERACK = '5'
SLOT_RESERVOIR = '3'
SLOT_TIPRACK_20 = '7'
SLOT_TIPRACK_300 = '11'          # 8-channel only, so partial-tip rules do not apply
SLOT_DILUTION_PLATE = '6'
# Slots 4, 8, 9 and 10 MUST BE EMPTY. Slot 9 in particular can never be a single-nozzle
# target: the fixed trash bin in slot 12 is directly north of it, and that is exactly
# what failed on 2026-09-04. See DECK RULES in W-O-SC-03.

ASSAY_PLATE_LOADNAME = 'nest_96_wellplate_200ul_flat'
DILUTION_PLATE_LOADNAME = 'nest_96_wellplate_200ul_flat'
RESERVOIR_LOADNAME = 'nest_12_reservoir_15ml'
TUBERACK_LOADNAME = '3d_printed_tuberack_1.5ml'
HS_ADAPTER_LOADNAME = 'opentrons_96_flat_bottom_adapter'

PLATE_CRASH_HEADSPACE_MM = 5.0
DISPENSE_ASPIRATE_CLEARANCE_MM = 13.0 + PLATE_CRASH_HEADSPACE_MM
SAMPLE_DISPENSE_CLEARANCE_MM = 1.0 + PLATE_CRASH_HEADSPACE_MM
DILUTION_DISPENSE_CLEARANCE_MM = 2.0
# Mixing depth. A CLEARANCE IS NOT A SUBMERSION DEPTH: 60 uL in a NEST 200 uL flat well
# stands only 1.63 mm deep, so mixing at the 2.0 mm dispense clearance draws pure air -
# which is exactly what this check found on 2026-09-04. Phase 5 now computes the depth
# the same way the real protocols do, so it re-tests that fix every time you run it.
DILUTION_MIX_MIN_Z_MM = 0.5
DILUTION_MIX_DEPTH_FRACTION = 0.5
DILUTION_BLOWOUT_ABOVE_MM = 1.0
DILUTION_FILL_VOL_UL = 60.0      # the real protocols' DILUTION_TOTAL_VOL_UL

# Droplet control, matching the real protocols. Blow out clears the residual inside the
# tip; touch tip sheds the droplet clinging to the outside. Watch for both working.
TOUCH_TIP_RADIUS = 0.8
TOUCH_TIP_V_OFFSET_MM = -1.5

WATER_WELL = 'A3'                # reservoir well holding the water
TEST_TUBE = 'A1'                 # tube rack well holding ~200 uL water

# Small volumes - this is a geometry check, not a dispensing-accuracy check.
SMALL_VOL_UL = 5.0               # p20 moves
P20_STROKE_UL = 20.0             # one full p20 stroke
COLUMN_VOL_UL = 100.0            # p300 8-channel move

# Heater-Shaker check. Short and hot enough to prove the setpoint is honoured;
# 40 degC is inside the module's 37-95 degC controllable range.
HS_TEMP_C = 40
HS_SLOW_RPM = 250                # the incubation speed
HS_FAST_RPM = 1000               # the development speed
HS_HOLD_MIN = 1


def run(protocol: protocol_api.ProtocolContext):

    def phase(number, title, detail):
        """Announce a phase and, if asked, wait for the operator before doing it."""
        protocol.comment('')
        protocol.comment('=' * 70)
        protocol.comment(f'  PHASE {number}: {title}')
        protocol.comment(f'  WATCH FOR: {detail}')
        protocol.comment('=' * 70)
        if STOP_BETWEEN_PHASES:
            protocol.pause(f'PHASE {number} - {title}. {detail} Resume when ready.')

    # =============================================================================
    # PHASE 1 - load everything. No movement, but this is where an adapter/plate
    # incompatibility or a missing custom labware definition fails, and it fails
    # instantly instead of mid-run.
    # =============================================================================
    protocol.comment('PHASE 1: loading labware and modules')
    hs_mod = protocol.load_module('heaterShakerModuleV1', HS_SLOT)
    # The adapter MUST be declared. Without it the API assumes the plate sits directly
    # on the module, undershoots every well's Z and drives the tip into the plate.
    hs_adapter = hs_mod.load_adapter(HS_ADAPTER_LOADNAME)
    assay_plate = hs_adapter.load_labware(ASSAY_PLATE_LOADNAME)
    dilution_plate = protocol.load_labware(DILUTION_PLATE_LOADNAME, SLOT_DILUTION_PLATE)
    reservoir = protocol.load_labware(RESERVOIR_LOADNAME, SLOT_RESERVOIR)
    tuberack = protocol.load_labware(TUBERACK_LOADNAME, SLOT_TUBERACK)
    tr_20 = protocol.load_labware('opentrons_96_tiprack_20ul', SLOT_TIPRACK_20)
    tr_300 = protocol.load_labware('opentrons_96_tiprack_300ul', SLOT_TIPRACK_300)

    p20 = protocol.load_instrument('p20_multi_gen2', 'right', tip_racks=[tr_20])
    p300 = protocol.load_instrument('p300_multi_gen2', 'left', tip_racks=[tr_300])

    water = reservoir[WATER_WELL]
    protocol.comment(f'  loaded OK. Water expected in reservoir {WATER_WELL} and '
                     f'tube {TEST_TUBE}.')

    hs_mod.close_labware_latch()
    protocol.comment('  latch closed')

    # Single-nozzle mode, anchored at H1 (the FRONT-most nozzle). With A1 the gantry
    # cannot reach the assay plate on the Heater-Shaker in slot 1.
    # Only the p20: the real protocols keep the p300 in 8-channel configuration for the
    # whole run, which is what allows its rack to sit in slot 11. Testing the p300 in
    # single-nozzle mode here would prove something the assay never does.
    p20.configure_nozzle_layout(style=SINGLE, start='H1', tip_racks=[tr_20])

    # =============================================================================
    # PHASE 2 - the move that crashed before: P20 single-nozzle into assay plate A1,
    # the lowest and most crash-prone target in the whole workflow.
    # =============================================================================
    if CHECK_P20_TO_PLATE:
        phase(2, 'P20 single-nozzle -> assay plate A1',
              f'tip should stop {SAMPLE_DISPENSE_CLEARANCE_MM:g} mm above the well '
              'bottom, not touch it. Then check NO DROPLET stays on the tip after the '
              'blow-out and touch-off.')
        p20.pick_up_tip()
        for well in ('A1', 'H1', 'A12', 'H12'):
            # The four corners: if the plate or adapter is off, a corner shows it worst.
            p20.aspirate(SMALL_VOL_UL, water)
            p20.touch_tip(water, radius=TOUCH_TIP_RADIUS,
                          v_offset=TOUCH_TIP_V_OFFSET_MM)
            p20.dispense(SMALL_VOL_UL,
                         assay_plate[well].bottom(z=SAMPLE_DISPENSE_CLEARANCE_MM))
            # Blow out + touch off, exactly as the real spike step does. Watch that no
            # droplet stays on the tip as it leaves the well.
            p20.blow_out(assay_plate[well].bottom(z=SAMPLE_DISPENSE_CLEARANCE_MM))
            p20.touch_tip(assay_plate[well], radius=TOUCH_TIP_RADIUS,
                          v_offset=TOUCH_TIP_V_OFFSET_MM)
            protocol.comment(f'  dispensed into {well}')
        p20.drop_tip()

    # =============================================================================
    # PHASE 3 - P20 reaching into a 1.5 mL tube on the custom 3D-printed rack. This is
    # the deepest descent of the run and depends entirely on the custom labware's
    # declared depth being right.
    # =============================================================================
    if CHECK_P20_FROM_TUBERACK:
        phase(3, 'P20 single-nozzle <- tube rack',
              'tip should enter the tube without hitting the bottom. This is the '
              'deepest move of the run.')
        p20.pick_up_tip()
        p20.aspirate(SMALL_VOL_UL, tuberack[TEST_TUBE])
        p20.dispense(SMALL_VOL_UL, assay_plate['A2'].bottom(z=SAMPLE_DISPENSE_CLEARANCE_MM))
        p20.drop_tip()

    # =============================================================================
    # PHASE 4 - P20 to/from the dilution plate in the back row.
    # =============================================================================
    if CHECK_P20_FROM_DILUTION:
        phase(4, 'P20 single-nozzle <-> dilution plate',
              'gantry should reach slot ' + SLOT_DILUTION_PLATE + ' without complaint.')
        p20.pick_up_tip()
        p20.aspirate(SMALL_VOL_UL, water)
        p20.dispense(SMALL_VOL_UL,
                     dilution_plate['A1'].bottom(z=DILUTION_DISPENSE_CLEARANCE_MM))
        p20.aspirate(SMALL_VOL_UL,
                     dilution_plate['A1'].bottom(z=DILUTION_DISPENSE_CLEARANCE_MM))
        p20.dispense(SMALL_VOL_UL, assay_plate['A3'].bottom(z=SAMPLE_DISPENSE_CLEARANCE_MM))
        p20.drop_tip()

    # =============================================================================
    # PHASE 5 - the p20 filling a dilution well the way U-O-00a does, then MIXING it.
    # This is the phase that caught the silent mixing bug: the tip must actually go
    # UNDER the liquid, and in a well this shallow that is a fraction of a millimetre.
    # Watch the tip and confirm liquid moves in and out of it.
    # =============================================================================
    if CHECK_P300_TO_DILUTION:
        well = dilution_plate['B1']
        area = math.pi * (well.diameter / 2) ** 2
        full_h = DILUTION_FILL_VOL_UL / area
        low_h = (DILUTION_FILL_VOL_UL - P20_STROKE_UL) / area
        mix_z = max(DILUTION_MIX_MIN_Z_MM, low_h * DILUTION_MIX_DEPTH_FRACTION)
        blow_z = min(full_h + DILUTION_BLOWOUT_ABOVE_MM, well.depth - 1.0)

        phase(5, 'P20 fill + MIX in the dilution plate  [the 2026-09-04 bug]',
              f'{DILUTION_FILL_VOL_UL:g} uL stands only {full_h:.2f} mm deep and drops '
              f'to {low_h:.2f} mm mid-stroke, so the tip mixes at {mix_z:.2f} mm. '
              'CONFIRM LIQUID ACTUALLY ENTERS THE TIP on every stroke - if it draws '
              'air, the real dilutions are not being mixed either.')
        protocol.comment(f'  liquid {full_h:.2f} mm deep, mixing tip at {mix_z:.2f} mm, '
                         f'blow-out at {blow_z:.2f} mm')
        p20.pick_up_tip()
        for _ in range(int(DILUTION_FILL_VOL_UL // P20_STROKE_UL)):
            p20.aspirate(P20_STROKE_UL, water)
            p20.dispense(P20_STROKE_UL,
                         well.bottom(z=DILUTION_DISPENSE_CLEARANCE_MM))
        p20.mix(3, P20_STROKE_UL, well.bottom(z=mix_z))
        p20.blow_out(well.bottom(z=blow_z))
        p20.drop_tip()

    # =============================================================================
    # PHASE 6 - back to full 8-channel for the whole-column reagent adds. Checks that
    # all eight nozzles clear the plate at the raised reagent clearance.
    # =============================================================================
    if CHECK_P300_MULTI_COLUMN:
        phase(6, 'P300 8-channel -> assay plate column 1',
              f'all 8 tips should stop {DISPENSE_ASPIRATE_CLEARANCE_MM:g} mm above the '
              'well bottoms, well clear of the liquid already there.')
        p300.configure_nozzle_layout(style=ALL, tip_racks=[tr_300])
        p300.pick_up_tip()
        for column in ('A1', 'A12'):
            p300.aspirate(COLUMN_VOL_UL, water)
            p300.dispense(COLUMN_VOL_UL,
                          assay_plate[column].bottom(z=DISPENSE_ASPIRATE_CLEARANCE_MM))
            p300.touch_tip()
            protocol.comment(f'  dispensed into column {column}')
        p300.drop_tip()

    # =============================================================================
    # PHASE 7 - the module itself: does it actually reach the setpoint, and does the
    # plate stay put at both shake speeds? A plate that walks at 1000 rpm here would
    # cross-contaminate the real assay.
    # =============================================================================
    if CHECK_HEATER_SHAKER:
        phase(7, 'Heater-Shaker',
              f'{HS_TEMP_C} degC then {HS_SLOW_RPM} and {HS_FAST_RPM} rpm. Watch that '
              'the plate does not shift in the latch and nothing splashes out.')
        hs_mod.set_and_wait_for_temperature(HS_TEMP_C)
        protocol.comment(f'  reached {HS_TEMP_C} degC')
        for rpm in (HS_SLOW_RPM, HS_FAST_RPM):
            hs_mod.set_and_wait_for_shake_speed(rpm)
            protocol.comment(f'  shaking at {rpm} rpm')
            protocol.delay(minutes=HS_HOLD_MIN)
        hs_mod.deactivate_shaker()
        hs_mod.deactivate_heater()
        protocol.comment('  shaker and heater off')

    # =============================================================================
    # PHASE 8 - release, so the plate can be lifted out.
    # =============================================================================
    hs_mod.open_labware_latch()
    protocol.comment('')
    protocol.comment('=' * 70)
    protocol.comment('  DECK CHECK COMPLETE - every move above succeeded.')
    protocol.comment('  Empty the plates and you are clear to run the real protocol.')
    protocol.comment('=' * 70)
