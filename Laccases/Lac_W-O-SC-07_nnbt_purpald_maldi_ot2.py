"""
W-O-SC-07  NNBT / Purpald laccase HTS + MALDI time course          Kristine Toft Johansen s215098
=======================================================================================
Dilutes every enzyme from its Bradford reading to one shared concentration, builds a
lactaldehyde standard curve, lays everything out in row-wise triplicates, runs the
assay and hands off to the plate reader at A530.

WELL RECIPE   10 uL spike + 140 uL reaction mix = 150 uL, 2 h @ 40 C,
              + 50 uL Purpald = 200 uL. That is the plate's full capacity.

STEPS
  U-O-00a  buffer into the dilution wells
  U-O-00b  enzyme stock in, mix            -> every enzyme at one concentration
  U-O-00c  build the lactaldehyde curve    -> serial or direct, script decides
  U-O-01   10 uL spikes into the assay plate, in plate order
  U-O-02   140 uL reaction mix, 8-channel, whole columns   -> reaction starts
  U-M-01   PAUSE - parafilm on
  U-O-03   incubate 40 C / 250 rpm / 120 min
  U-M-02   PAUSE - parafilm off
  U-O-04   50 uL Purpald, 8-channel
  U-O-05   develop 40 C / 1000 rpm / 10 min
  U-C-06   latch open, read A530

TO RUN A BATCH   edit SECTION 1 only, then:
  python3 <this file>                      preview: plate map, recipes, what to load
  opentrons_simulate -L Labware <this file>   catches most errors, NOT deck slot 9
  upload to the Opentrons app

TWO RULES THAT HAVE ALREADY COST A RUN
  1  Slot 9 can never take a single-nozzle move: the trash in slot 12 is north of it,
     and opentrons_simulate does not model that. check_deck() below enforces it.
  2  A HEIGHT is not a DEPTH. 60 uL in these wells stands 1.63 mm deep, so a tip at
     2 mm draws air. Anything that must reach liquid is COMPUTED, never guessed.
"""

import math

# The opentrons package only exists on the robot / in the simulator. Guarding the
# import lets you run this same file locally (python3 <file>) to preview the layout
# and the volume budget before you ever touch the robot - see the __main__ block.
try:
    from opentrons import protocol_api
    from opentrons.protocol_api import SINGLE, ALL
except ImportError:                                  # pragma: no cover - local preview
    protocol_api = None
    SINGLE = ALL = None

def add_parameters(parameters):
    """Shown in the Opentrons app at run setup, so a part-used rack or a MALDI run
    needs no edit to this file."""
    parameters.add_str(
        display_name='20 uL tip row', variable_name='tip20_row', default='A',
        choices=[{'display_name': r, 'value': r} for r in 'ABCDEFGH'],
        description='Row of the first UNUSED tip in the 20 uL rack.')
    parameters.add_int(
        display_name='20 uL tip column', variable_name='tip20_col', default=1,
        minimum=1, maximum=12,
        description='Column of the first UNUSED tip in the 20 uL rack.')
    parameters.add_int(
        display_name='300 uL tip column', variable_name='tip300_col', default=1,
        minimum=1, maximum=12,
        description='First UNUSED column in the 300 uL rack (8-channel takes columns).')
    parameters.add_bool(
        display_name='MALDI time course', variable_name='maldi_enabled', default=False,
        description='Pause the incubation at intervals and spot the MALDI target.')
    parameters.add_int(
        display_name='MALDI interval (min)', variable_name='maldi_interval', default=30,
        minimum=5, maximum=120,
        description='Minutes between spotting rounds. Ignored unless MALDI is on.')


metadata = {
    'apiLevel': '2.20',
    'protocolName': 'W-O-SC-07 - NNBT Purpald laccase HTS (auto-dilution + standard curve)',
    'description': (
        'Dilutes each enzyme from its Bradford concentration to a common assay '
        'concentration, builds a lactaldehyde standard curve, lays enzymes and '
        'controls out in triplicate on one 96-well plate, adds reaction mix, '
        'incubates 2 h at 40 degC, adds Purpald, develops and hands off for A530.'
    ),
}


# =======================================================================================
# SECTION 1 - THE BATCH   *** edit this and nothing else between runs ***
# =======================================================================================
# Tubes go on the rack in THIS ORDER: A1, B1, C1, D1, A2, ... The list length is the
# batch size. 'conc' is your Bradford number in any unit, as long as the target below
# uses the same one - only ratios are ever taken.
ENZYME_BATCH = [
    {'name': 'Lac-01', 'conc': 2.40},          # tube A1
    {'name': 'Lac-02', 'conc': 1.85},          # tube B1
    {'name': 'Lac-03', 'conc': 3.10},          # tube C1
]

ENZ_TARGET_CONC = None                 # concentration every enzyme is diluted to; None = auto
ENZ_TARGET_AUTO_FRACTION = 0.90        # auto = this x the weakest enzyme, so all get diluted

HEAT_INACT_NAME = 'Heat-inactivated'   # negative control: killed enzyme, tube C6
HEAT_INACT_CONC = 2.40                 # its own Bradford reading; heat can precipitate protein
INCLUDE_HEAT_INACT = True              # False + ENZYME_BATCH=[] gives a curve-only pilot plate

# --- lactaldehyde positive control -----------------------------------------------------
LAC_STOCK_UM = 100000.0                # tube D6. 1 M diluted 1:10 = 100000 uM, NOT 10000
LAC_GRADIENT_UM = [500, 1500, 2500, 3500, 4500, 6000]   # uM IN THE 150 uL REACTION (x-axis)
GRADIENT_STRATEGY = 'serial'           # 'serial' = biggest, most accurate volumes, errors
                                       # compound; 'direct' = off the stock, independent
                                       # errors, smaller volumes. Script flags each step.
DILUTION_WELL_VOL_UL = 100.0           # volume made up in EACH dilution well (3 x 10 uL used)


# =======================================================================================
# SECTION 2 - ASSAY CHEMISTRY   (fixed unless the assay itself changes)
# =======================================================================================
SPIKE_VOL_UL = 10.0                    # sample into each well; below the p300's 20 uL floor
REACTION_MIX_VOL_UL = 140.0            # buffer + NNBT, premixed off-deck
PURPALD_VOL_UL = 50.0                  # 10 + 140 + 50 = 200 uL = the plate's capacity
REACTION_VOL_UL = SPIKE_VOL_UL + REACTION_MIX_VOL_UL    # 150 uL during the incubation

REPLICATES = 3                         # triplicates
REPLICATE_ORIENTATION = 'row'          # 'row' = A1,A2,A3 side by side; 'column' = A1,B1,C1
WELLS_PER_COLUMN = 8
ROW_LETTERS = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
PLATE_WELL_COUNT = 96

INCUBATION_TEMP_C = 40                 # inside the module's 37-95 C controllable range
INCUBATION_RPM = 250
INCUBATION_MIN = 120

DEVELOP_TEMP_C = 40
DEVELOP_RPM = 1000                     # not the 3000 ceiling: parafilm is off, wells brim-full
DEVELOP_MIN = 10

REACTION_MIX_FLOW_RATE_SCALE = 0.5     # NNBT is in acetonitrile and drips at full speed

P20_FLOW_RATE_UL_S = 15.0              # default 7.6, hardware max 24
P300_FLOW_RATE_UL_S = 150.0            # default 94, hardware max 275

# =======================================================================================
# SECTION 3 - HARDWARE
# =======================================================================================
# DECK - REARRANGED FROM W-O-SC-05. Both pipettes now work single-nozzle, so BOTH tip
# racks must be reachable in partial-tip mode, which W-O-SC-05's deck could not do.
# The rule: in single-nozzle mode the pipette body overhangs the slot to the NORTH
# (slot + 3), so nothing tall may sit there. Tall = tip racks (64 mm), tube rack (54 mm);
# a plate (14 mm) or the reservoir (31 mm) is fine. check_deck() enforces it.
HS_SLOT = 1                            # Heater-Shaker + assay plate  (north: 4, empty)
SLOT_TUBERACK = '3'                    # enzyme + lactaldehyde stocks (north: 6, plate 14mm)
SLOT_TIPRACK_20 = '5'                  # 20 uL tips                   (north: 8, reservoir)
SLOT_DILUTION_PLATE = '6'              # every dilution is made here  (north: 9, empty)
SLOT_TIPRACK_300 = '7'                 # 300 uL tips                  (north: 10, empty)
SLOT_RESERVOIR = '8'                   # reaction mix / buffer / Purpald (north: 11, empty)
# Slots 2, 4, 9, 10, 11 MUST STAY EMPTY - they are clearance, not spare space.
# With MALDI on, the target takes over slot 3 once the tube rack's job is done: there is
# no free slot a single-nozzle move can reach, and the two are never needed at once.
# Slot 9 can never hold single-nozzle work at all: the trash in slot 12 is north of it,
# and opentrons_simulate does not model that bin. Slot 2 is refused as Heater-Shaker
# adjacent. THE ONE UNPROVEN NEIGHBOUR here is the 31 mm reservoir north of the 20 uL
# tip rack - it simulates clean, but watch that pickup on the first real run.

# PARTLY USED TIP RACKS. Set to the first UNUSED tip and the run starts there instead of
# A1. For the 8-channel p300 give a top-row well ('A3') - it takes whole columns.
# Tip start positions are RUNTIME PARAMETERS - set them in the Opentrons app at run
# setup instead of editing this file. The values below are only the app's defaults.
FIRST_TIP_20 = 'A1'                    # first unused 20 uL tip
FIRST_TIP_300 = 'A1'                   # first unused 300 uL column (8-ch takes whole columns)

ASSAY_PLATE_LOADNAME = 'nest_96_wellplate_200ul_flat'
DILUTION_PLATE_LOADNAME = 'nest_96_wellplate_200ul_flat'
RESERVOIR_LOADNAME = 'nest_12_reservoir_15ml'
TUBERACK_LOADNAME = '3d_printed_tuberack_1.5ml'
HS_ADAPTER_LOADNAME = 'opentrons_96_flat_bottom_adapter'   # MUST be declared or Z is wrong

TUBERACK_ROWS_PER_COLUMN = 4           # rack is 4 rows (A-D) x 6 columns = 24 tubes
TUBERACK_COLUMNS = 6
TUBERACK_CAPACITY = TUBERACK_ROWS_PER_COLUMN * TUBERACK_COLUMNS
TUBE_HEAT_INACT = 'C6'                 # fixed, so it never moves between runs
TUBE_LACTALDEHYDE = 'D6'               # fixed
MAX_ENZYME_TUBES = TUBERACK_CAPACITY - 2                # the real batch-size ceiling

RESERVOIR_REACTION_MIX_WELLS = ['A1', 'A2']   # two wells: a full plate needs >13 mL
RESERVOIR_BUFFER_WELL = 'A3'                  # dilution buffer + the blank's spike
RESERVOIR_PURPALD_WELLS = ['A4']
RESERVOIR_USABLE_VOL_UL = 13_000.0     # 15 mL nominal, minus fill margin
RESERVOIR_DEAD_VOL_UL = 1_000.0        # pour this much extra so tips never hit air

P20_MIN_UL, P20_MAX_UL = 1.0, 20.0     # p20 working range
P300_MIN_UL, P300_MAX_UL = 20.0, 300.0 # p300 working range; also does single >20 uL

# --- tip heights -----------------------------------------------------------------------
# HEIGHTS are absolute, above the well bottom, and are used where the tip must stay OUT
# of the liquid. Anything that must reach INTO liquid is computed at run time instead
# (see mix_location / aspirate_depth), because 100 uL in these wells is only 2.7 mm deep.
EXTRA_HEIGHT_MARGIN_MM = 2.0           # padding for the residual real-vs-modelled Z gap
REAGENT_DISPENSE_HEIGHT_MM = 13.0 + EXTRA_HEIGHT_MARGIN_MM   # 8-channel reagents onto liquid
SPIKE_DISPENSE_HEIGHT_MM = 1.0 + EXTRA_HEIGHT_MARGIN_MM      # spikes into empty assay wells
DILUTION_DISPENSE_HEIGHT_MM = 0.5      # into a dilution well; must stay under 2.7 mm

DILUTION_MIX_REPS = 3                  # 40 uL x 3 = 120 uL total mix volume
DILUTION_MIX_STROKE_FRACTION = 0.4     # of the well volume; smaller stroke keeps the
                                       # surface high enough for the tip to stay under it
DILUTION_MIX_MIN_HEIGHT_MM = 1.0       # floor: a p300 tip SEALS on a flat well bottom
                                       # below ~0.8 mm and then moves no liquid at all
DILUTION_MIX_DEPTH_FRACTION = 0.5      # of the liquid left at the BOTTOM of a mix stroke
DILUTION_ASPIRATE_MIN_HEIGHT_MM = 0.8  # same floor problem, same fix
DILUTION_ASPIRATE_DEPTH_FRACTION = 0.5 # of the liquid left AFTER the aspirate stroke
BLOWOUT_ABOVE_SURFACE_MM = 1.0         # blow out above the surface - under it foams protein

# --- MALDI time course -----------------------------------------------------------------
# Switched on, and the interval set, in the Opentrons app. During the incubation the
# robot stops the shaker every interval, draws a few uL from well 1 of each condition
# and spots it in triplicate. Wells 2 and 3 are never touched.
#
# *** REPLACE MALDI_LOADNAME AND THE GRID BELOW with your colleague's definition. ***
# The default is a placeholder 384-spot target (16 x 24, 4.5 mm pitch) whose Z is a
# guess - it will run, but do not trust its heights until you have measured yours.
# A bare MALDI target is not SBS footprint; it needs a holder to sit in a deck slot.
MALDI_LOADNAME = 'maldi_384_target_sbs_holder'
MALDI_ROWS = list('ABCDEFGHIJKLMNOP')  # 16 rows
MALDI_COLUMNS = 24
MALDI_SLOT = SLOT_TUBERACK             # the target REPLACES the tube rack (see below)
MALDI_SPOT_VOL_UL = 1.0                # p20 minimum; MALDI wants a small dry spot
MALDI_SPOT_REPLICATES = 3              # spotting replicates from ONE aspirate, not
                                       # reaction replicates - those stay on the plate
MALDI_ASPIRATE_EXTRA_UL = 1.0          # surplus so the last spot is not the piston's end
MALDI_SPOT_HEIGHT_MM = 1.0             # above the spot surface - TUNE THIS on hardware
MALDI_SAMPLE_HEIGHT_MM = 3.0           # draw from here in the 150 uL assay well

# --- droplet control -------------------------------------------------------------------
# A hanging drop is both a short delivery and a carry-over route. Blow-out clears the
# INSIDE of the tip; touch-tip sheds the drop on the OUTSIDE, against the source after
# aspirating and the destination after dispensing.
BLOW_OUT_ENABLED = True
TOUCH_TIP_ENABLED = False
TOUCH_TIP_RADIUS = 0.8                 # fraction of well radius; <1 keeps clear of the wall
TOUCH_TIP_V_OFFSET_MM = -1           # below the well top, so above the liquid everywhere


# =======================================================================================
# SECTION 4 - LAYOUT ENGINE   (derived; you should not need to edit any of it)
# =======================================================================================
REGION_BUFFER_CTRL = 'neg_buffer'
REGION_HEAT_CTRL = 'neg_heat'
REGION_POSITIVE = 'pos_lactaldehyde'
REGION_SAMPLE = 'enzyme'

# Colours used for the Opentrons app's plate map and the local HTML preview.
REGION_COLORS = {
    REGION_BUFFER_CTRL: '#78706f',   # warm grey - the blank
    REGION_HEAT_CTRL: '#3f6d8f',     # steel blue - the killed-enzyme control
    REGION_POSITIVE: '#c2621a',      # amber - the standard curve
    REGION_SAMPLE: '#2f7d5f',        # green - the enzymes under test
}
REGION_LABELS = {
    REGION_BUFFER_CTRL: 'NEG - buffer only',
    REGION_HEAT_CTRL: 'NEG - heat-inactivated enzyme',
    REGION_POSITIVE: 'POS - lactaldehyde standard',
    REGION_SAMPLE: 'Enzyme sample',
}


# --- deck preflight -------------------------------------------------------------------
# Heights in mm, from the labware definitions. Anything at or above ~50 mm reliably
# blocks a single-nozzle move from the slot to its south; the fixed trash bin is taller
# still. Kept here so the check below is self-contained and readable.
LABWARE_HEIGHT_MM = {
    'opentrons_96_tiprack_20ul': 64.7,
    'opentrons_96_tiprack_300ul': 64.5,
    '3d_printed_tuberack_1.5ml': 54.0,
    'nest_12_reservoir_15ml': 31.4,
    'nest_96_wellplate_200ul_flat': 14.3,
    'TRASH': 999.0,
}
TRASH_SLOT = '12'
HS_ADJACENT_SLOTS = ('2',)              # beside the Heater-Shaker in slot 1
OUT_OF_REACH_SLOTS = ('10', '11', '12')  # unreachable with the 'H1' anchor
BLOCKING_HEIGHT_MM = 50.0


def check_mix_geometry():
    """Refuse to load if a computed tip height is outside the liquid.

    Both failure modes have already happened: a tip parked ABOVE the surface draws air
    and silently leaves the dilutions unmixed; a tip parked ON a flat well bottom seals
    the orifice and moves nothing. Both look fine in the run log. This catches either
    before the robot moves.
    """
    area = math.pi * (6.85 / 2) ** 2                 # NEST 200 uL flat well
    V = DILUTION_WELL_VOL_UL
    stroke = min(V * DILUTION_MIX_STROKE_FRACTION, P300_MAX_UL)
    surface = (V - stroke) / area                    # level at the bottom of a stroke
    z = max(DILUTION_MIX_MIN_HEIGHT_MM, surface * DILUTION_MIX_DEPTH_FRACTION)
    if z >= surface:
        raise ValueError(
            f'mixing tip would sit at {z:.2f} mm but a {stroke:g} uL stroke drops the '
            f'surface to {surface:.2f} mm - it would draw AIR. Lower '
            'DILUTION_MIX_MIN_HEIGHT_MM or DILUTION_MIX_STROKE_FRACTION.')
    if z < 0.8:
        raise ValueError(
            f'mixing tip would sit at {z:.2f} mm above a FLAT well bottom - a p300 tip '
            'seals below ~0.8 mm and moves no liquid. Raise DILUTION_MIX_MIN_HEIGHT_MM.')


check_mix_geometry()


def check_deck():
    """Refuse to load on a deck rule we have already been burned by.

    opentrons_simulate passed the layout the robot then refused (dilution plate in slot
    9, trash north of it) because it does not model the bin. This does, so a bad slot
    fails on your laptop instead of 40 s into a run. Not a substitute for simulating,
    and neither replaces running W-O-SC-90 on the hardware."""
    occupied = {
        str(HS_SLOT): ASSAY_PLATE_LOADNAME,
        SLOT_TUBERACK: TUBERACK_LOADNAME,
        SLOT_RESERVOIR: RESERVOIR_LOADNAME,
        SLOT_DILUTION_PLATE: DILUTION_PLATE_LOADNAME,
        SLOT_TIPRACK_20: 'opentrons_96_tiprack_20ul',
        SLOT_TIPRACK_300: 'opentrons_96_tiprack_300ul',
        TRASH_SLOT: 'TRASH',
    }
    if len(occupied) < 7:
        raise ValueError('two pieces of labware are assigned to the same deck slot.')

    # Slots the pipette enters in single-nozzle mode. The p300 is 8-channel throughout,
    # so its rack is deliberately NOT in this list.
    single_nozzle_targets = {
        str(HS_SLOT): 'assay plate', SLOT_TUBERACK: 'tube rack',
        SLOT_RESERVOIR: 'reservoir', SLOT_DILUTION_PLATE: 'dilution plate',
        SLOT_TIPRACK_20: '20 uL tips',
        SLOT_TIPRACK_300: '300 uL tips',
    }
    problems = []
    for slot, what in single_nozzle_targets.items():
        if slot in OUT_OF_REACH_SLOTS:
            problems.append(f'{what} is in slot {slot}, which the H1 nozzle anchor '
                            'cannot reach at all. Use slots 1-9.')
        if slot in HS_ADJACENT_SLOTS:
            problems.append(f'{what} is in slot {slot}, directly beside the '
                            'Heater-Shaker - refused for any multi-channel pipette.')
        north = str(int(slot) + 3)
        if north == TRASH_SLOT:
            problems.append(
                f'{what} is in slot {slot}, and the fixed trash bin in slot 12 is '
                'directly north of it. This is the exact failure that killed the run '
                'of 2026-09-04. opentrons_simulate will NOT catch this. Move it.')
        elif north in occupied:
            h = LABWARE_HEIGHT_MM.get(occupied[north], 0.0)
            if h >= BLOCKING_HEIGHT_MM:
                problems.append(
                    f'{what} is in slot {slot} with {occupied[north]} ({h:g} mm) north '
                    f'of it in slot {north}. Anything this tall blocks the move.')
    if problems:
        raise ValueError('DECK RULE VIOLATION (see DECK RULES in SECTION 3):\n  - '
                         + '\n  - '.join(problems))


check_deck()


PLATE_COLUMN_COUNT = PLATE_WELL_COUNT // WELLS_PER_COLUMN     # 12


def triplicate_groups(count):
    """The first `count` triplicate groups, as lists of well names.

    'row': replicates side by side (A1,A2,A3), handed out one 3-column BLOCK at a time,
    top to bottom, before moving right. Filling row A across first would leave every
    column part-empty and make a handful of samples cost a full 12-column reagent add.
    Side effect: columns 1-3 end up all controls and standards.
    'column': replicates stacked (A1,B1,C1), packed continuously.
    """
    if REPLICATE_ORIENTATION == 'row':
        if PLATE_COLUMN_COUNT % REPLICATES:
            raise ValueError(
                f"REPLICATE_ORIENTATION 'row' needs REPLICATES ({REPLICATES}) to divide "
                f'the plate\'s {PLATE_COLUMN_COUNT} columns evenly. Use '
                "REPLICATE_ORIENTATION = 'column' instead.")
        groups = [[f'{row}{block * REPLICATES + r + 1}' for r in range(REPLICATES)]
                  for block in range(PLATE_COLUMN_COUNT // REPLICATES)
                  for row in ROW_LETTERS]
    elif REPLICATE_ORIENTATION == 'column':
        names = [f'{row}{col + 1}'
                 for col in range(PLATE_COLUMN_COUNT) for row in ROW_LETTERS]
        groups = [names[i:i + REPLICATES] for i in range(0, len(names), REPLICATES)]
    else:
        raise ValueError("REPLICATE_ORIENTATION must be 'row' or 'column'")
    return groups[:count]


def tube_well_name(index):
    """Tube-rack well for the `index`-th enzyme (0-based): 4 rows per column, A-D."""
    return f'{ROW_LETTERS[index % TUBERACK_ROWS_PER_COLUMN]}{index // TUBERACK_ROWS_PER_COLUMN + 1}'


def dilution_well_name(index):
    """Dilution-plate well for the `index`-th dilution (0-based), column-then-row."""
    return f'{ROW_LETTERS[index % WELLS_PER_COLUMN]}{index // WELLS_PER_COLUMN + 1}'


def resolve_target_conc():
    """The common concentration every enzyme is diluted down to."""
    if ENZ_TARGET_CONC is not None:
        return float(ENZ_TARGET_CONC)
    # In pilot mode there may be no enzymes at all; the heat-inactivated control is then
    # the only protein on the deck, so it sets the target. With neither, nothing is
    # diluted and the value is unused.
    measured = [e['conc'] for e in ENZYME_BATCH]
    if INCLUDE_HEAT_INACT:
        measured.append(HEAT_INACT_CONC)
    if not measured:
        return 0.0
    return min(measured) * ENZ_TARGET_AUTO_FRACTION


def dilution_recipe(stock_conc, target_conc, total_vol=None):
    """(stock_uL, buffer_uL) making `total_vol` at `target_conc`. Raises rather than
    silently producing a dilution the pipettes cannot actually deliver."""
    total_vol = DILUTION_WELL_VOL_UL if total_vol is None else total_vol
    if stock_conc <= 0:
        raise ValueError('stock concentration must be > 0')
    stock_vol = round(total_vol * target_conc / stock_conc, 1)
    if stock_vol > total_vol + 1e-9:
        raise ValueError(
            f'stock at {stock_conc:g} is WEAKER than the target {target_conc:g}: it '
            f'would need {stock_vol:.1f} uL of stock in a {total_vol:g} uL well. Lower '
            'ENZ_TARGET_CONC (or ENZ_TARGET_AUTO_FRACTION), or drop that sample.'
        )
    if stock_vol < P20_MIN_UL:
        raise ValueError(
            f'stock at {stock_conc:g} is so concentrated that only {stock_vol:.2f} uL '
            f'would be needed - below the p20 minimum of {P20_MIN_UL:g} uL. Either '
            'raise DILUTION_WELL_VOL_UL, or pre-dilute that stock by hand and enter '
            'the pre-diluted concentration.'
        )
    return stock_vol, round(total_vol - stock_vol, 1)


def build_gradient_plan(dilution_index_start):
    """Plan the lactaldehyde curve, most concentrated point first so a point is always
    finished before anything dilutes from it. Each point is a dilution well holding the
    concentration that gives the requested in-reaction value once the 10 uL spike is
    diluted into 140 uL. GRADIENT_STRATEGY picks the source; only sources needing at
    least P20_MIN_UL of themselves are considered."""
    if GRADIENT_STRATEGY not in ('serial', 'direct'):
        raise ValueError("GRADIENT_STRATEGY must be 'serial' or 'direct'")
    spike_factor = REACTION_VOL_UL / SPIKE_VOL_UL     # 150/10 = 15x
    points = sorted(LAC_GRADIENT_UM, reverse=True)
    # (kind, address, concentration-in-that-vessel)
    sources = [('tube', TUBE_LACTALDEHYDE, float(LAC_STOCK_UM))]
    plan = []
    for i, final_um in enumerate(points):
        spike_um = final_um * spike_factor
        usable = [s for s in sources
                  if DILUTION_WELL_VOL_UL * spike_um / s[2] >= P20_MIN_UL
                  and s[2] >= spike_um]
        if not usable:
            # Two very different causes, so say which one it is. Reporting the wrong
            # one sends you looking at volumes when the real problem is the stock.
            strongest = max(s[2] for s in sources)
            if strongest < spike_um:
                raise ValueError(
                    f'gradient point {final_um:g} uM is STRONGER THAN THE STOCK CAN '
                    f'MAKE. It needs a {spike_um:g} uM spike (x{spike_factor:g} the '
                    f'in-reaction concentration, because 10 uL is diluted into 150 uL), '
                    f'but the most concentrated source available is {strongest:g} uM.\n'
                    f'    With LAC_STOCK_UM = {LAC_STOCK_UM:g} the highest reachable '
                    f'point is {LAC_STOCK_UM / spike_factor:.0f} uM in the reaction.\n'
                    f'    To reach {final_um:g} uM you need a stock of at least '
                    f'{spike_um:g} uM. Check the units: 1 M = 1e6 uM, so a 1:10 '
                    'dilution of a 1 M stock is 100000 uM, not 10000.')
            raise ValueError(
                f'gradient point {final_um:g} uM cannot be reached: every source '
                f'concentrated enough would need less than the p20 minimum of '
                f'{P20_MIN_UL:g} uL. Add a higher point to serial-dilute from, or '
                'raise DILUTION_WELL_VOL_UL.')
        # 'serial': least concentrated usable source = biggest, most accurate transfer.
        # 'direct': most concentrated usable source = the stock tube whenever it is
        #           pipettable, so no point inherits the previous point's error.
        source = (min(usable, key=lambda s: s[2]) if GRADIENT_STRATEGY == 'serial'
                  else max(usable, key=lambda s: s[2]))
        src_vol, buf_vol = dilution_recipe(source[2], spike_um)
        well = dilution_well_name(dilution_index_start + i)
        plan.append({
            'name': f'Lac-std {final_um:g} uM',
            'final_um': final_um,
            'spike_um': spike_um,
            'dilution_well': well,
            'source_kind': source[0],
            'source_address': source[1],
            'source_conc': source[2],
            'source_vol': src_vol,
            'buffer_vol': buf_vol,
            'serial': source[0] == 'dilution',
        })
        sources.append(('dilution', well, spike_um))
    return plan


def assign_reservoir_wells(vol_per_column, columns, well_names, reagent):
    """Which reservoir well each 8-channel column add draws from, plus how much to pour
    into each. A full plate needs 13.44 mL of reaction mix, more than one 15 mL well, so
    the load is spread EVENLY over the fewest wells that hold it. The setup sheet and the
    run both read this, so they cannot disagree."""
    total = vol_per_column * columns
    wells_needed = max(1, math.ceil(total / RESERVOIR_USABLE_VOL_UL))
    if wells_needed > len(well_names):
        raise ValueError(
            f'{reagent} needs {total / 1000:.1f} mL, i.e. {wells_needed} reservoir '
            f'wells at {RESERVOIR_USABLE_VOL_UL / 1000:g} mL each, but only '
            f'{len(well_names)} are configured. Add another well name to its '
            'RESERVOIR_*_WELLS list.')
    columns_per_well = math.ceil(columns / wells_needed)
    per_column, needs = [], {}
    for c in range(columns):
        well = well_names[c // columns_per_well]
        needs[well] = needs.get(well, 0.0) + vol_per_column
        per_column.append(well)
    return per_column, needs


def build_layout():
    """The single source of truth: dilution plan, plate assignment, volume budget.
    Pure Python, no opentrons objects, so it runs on your laptop."""
    n_enz = len(ENZYME_BATCH)
    if n_enz > MAX_ENZYME_TUBES:
        raise ValueError(
            f'{n_enz} enzymes but the tube rack only has {MAX_ENZYME_TUBES} free slots '
            f'({TUBERACK_CAPACITY} tubes minus the heat-inactivated and lactaldehyde '
            'stock tubes).')

    target = resolve_target_conc()

    # ---- dilution plate: enzymes, then heat-inactivated, then the gradient ----------
    dilutions = []
    for i, enz in enumerate(ENZYME_BATCH):
        stock_vol, buf_vol = dilution_recipe(enz['conc'], target)
        dilutions.append({
            'name': enz['name'],
            'region': REGION_SAMPLE,
            'dilution_well': dilution_well_name(i),
            'source_kind': 'tube',
            'source_address': tube_well_name(i),
            'source_conc': enz['conc'],
            'source_vol': stock_vol,
            'buffer_vol': buf_vol,
            'final_conc': target,
        })

    heat_dilution = None
    if INCLUDE_HEAT_INACT:
        hi_stock, hi_buf = dilution_recipe(HEAT_INACT_CONC, target)
        heat_dilution = {
            'name': HEAT_INACT_NAME,
            'region': REGION_HEAT_CTRL,
            'dilution_well': dilution_well_name(n_enz),
            'source_kind': 'tube',
            'source_address': TUBE_HEAT_INACT,
            'source_conc': HEAT_INACT_CONC,
            'source_vol': hi_stock,
            'buffer_vol': hi_buf,
            'final_conc': target,
        }

    gradient = build_gradient_plan(n_enz + (1 if INCLUDE_HEAT_INACT else 0))
    for g in gradient:
        g['region'] = REGION_POSITIVE

    # ---- plate assignment: controls first (fixed addresses), then samples ----------
    # entry = (region, label, source-of-the-10-uL-spike)
    entries = [(REGION_BUFFER_CTRL, 'Buffer only (0 uM)',
                ('reservoir', RESERVOIR_BUFFER_WELL))]
    if heat_dilution:
        entries.append((REGION_HEAT_CTRL, HEAT_INACT_NAME,
                        ('dilution', heat_dilution['dilution_well'])))
    entries += [(REGION_POSITIVE, g['name'], ('dilution', g['dilution_well']))
                for g in gradient]
    entries += [(REGION_SAMPLE, d['name'], ('dilution', d['dilution_well']))
                for d in dilutions]

    groups_needed = len(entries)
    max_groups = PLATE_WELL_COUNT // REPLICATES          # 32
    if groups_needed > max_groups:
        raise ValueError(
            f'{groups_needed} triplicate groups requested but only {max_groups} fit on '
            f'one plate. Room for {max_groups - (groups_needed - n_enz)} enzymes with '
            f'{len(gradient)} gradient points and '
            f'{2 if heat_dilution else 1} negative control(s).')

    groups = triplicate_groups(groups_needed)
    wells = [w for g in groups for w in g]
    assignments = []
    plate_map = {}
    for i, (region, label, source) in enumerate(entries):
        group = groups[i]
        assignments.append({'region': region, 'label': label,
                            'source': source, 'wells': group})
        for r, w in enumerate(group):
            plate_map[w] = {'region': region, 'label': label, 'replicate': r + 1}

    # The rightmost column holding a sample. Every column up to it is fully or partly
    # occupied, and the 8-channel reagent adds cover columns 1..columns_used.
    columns_used = max(int(w[1:]) for w in wells)

    # ---- volume budget -------------------------------------------------------------
    # Reagent adds are whole-column and 8-channel, so a partly-filled last column still
    # costs a full column's worth of reagent. assign_reservoir_wells() decides which
    # reservoir well each column is drawn from, and the SAME list drives both the
    # "what to load" sheet and the run itself - so what you are told to pour is exactly
    # what gets aspirated.
    rxn_wells, rxn_needs = assign_reservoir_wells(
        REACTION_MIX_VOL_UL * WELLS_PER_COLUMN, columns_used,
        RESERVOIR_REACTION_MIX_WELLS, 'reaction mix')
    purpald_wells, purpald_needs = assign_reservoir_wells(
        PURPALD_VOL_UL * WELLS_PER_COLUMN, columns_used,
        RESERVOIR_PURPALD_WELLS, 'Purpald')
    rxn_total = REACTION_MIX_VOL_UL * WELLS_PER_COLUMN * columns_used
    purpald_total = PURPALD_VOL_UL * WELLS_PER_COLUMN * columns_used
    all_dilutions = dilutions + ([heat_dilution] if heat_dilution else [])
    buffer_total = (sum(d['buffer_vol'] for d in all_dilutions)
                    + sum(g['buffer_vol'] for g in gradient)
                    + SPIKE_VOL_UL * REPLICATES)
    lac_stock_needed = sum(g['source_vol'] for g in gradient if not g['serial'])

    return {
        'target_conc': target,
        'dilutions': dilutions,
        'heat_dilution': heat_dilution,
        'all_dilutions': all_dilutions,
        'gradient': gradient,
        'assignments': assignments,
        'plate_map': plate_map,
        'wells_used': len(wells),
        'columns_used': columns_used,
        'plate_columns': [f'A{i + 1}' for i in range(columns_used)],
        'free_enzyme_slots': max_groups - groups_needed,
        'rxn_well_per_column': rxn_wells,
        'purpald_well_per_column': purpald_wells,
        'reservoir_needs': {**rxn_needs, **purpald_needs,
                            RESERVOIR_BUFFER_WELL: buffer_total},
        'volumes': {
            'reaction_mix': rxn_total,
            'purpald': purpald_total,
            'buffer': buffer_total,
            'lactaldehyde_stock': lac_stock_needed,
        },
    }


LAYOUT = build_layout()


# =======================================================================================
# SECTION 5 - TEXT RENDERING  (used both in the run log and in the local preview)
# =======================================================================================
REGION_SYMBOL = {
    REGION_BUFFER_CTRL: 'B',
    REGION_HEAT_CTRL: 'H',
    REGION_POSITIVE: 'P',
    REGION_SAMPLE: 'E',
}


def render_plate_map(layout):
    """ASCII plate map: one cell per well, symbol + group index."""
    index_of = {}
    for i, a in enumerate(layout['assignments']):
        for w in a['wells']:
            index_of[w] = i + 1
    lines = ['     ' + ''.join(f'{c + 1:>6}' for c in range(12))]
    for row in ROW_LETTERS:
        cells = []
        for col in range(1, 13):
            w = f'{row}{col}'
            info = layout['plate_map'].get(w)
            cells.append(f'{REGION_SYMBOL[info["region"]]}{index_of[w]:02d}'.rjust(6)
                         if info else '     .')
        lines.append(f'  {row}  ' + ''.join(cells))
    lines.append('')
    lines.append('  legend:  ' + '   '.join(
        f'{REGION_SYMBOL[r]}={REGION_LABELS[r]}' for r in REGION_SYMBOL))
    lines.append('  group numbers below; "." = empty well')
    return lines


def render_group_key(layout):
    lines = ['  #   wells            content']
    for i, a in enumerate(layout['assignments']):
        lines.append(f'  {i + 1:02d}  {",".join(a["wells"]):<16} {a["label"]}')
    return lines


def render_dilutions(layout):
    lines = []
    if layout['all_dilutions']:
        lines += [f'  every enzyme is diluted to {layout["target_conc"]:.3g} '
                  f'(same unit as the Bradford input), {DILUTION_WELL_VOL_UL:g} uL '
                  'per well',
                  '  dil.well  from        stock uL  buffer uL  content']
    else:
        lines.append('  no enzymes in this batch - standard curve only (pilot mode)')
    for d in layout['all_dilutions']:
        lines.append(f'  {d["dilution_well"]:<9} tube {d["source_address"]:<7} '
                     f'{d["source_vol"]:>8.1f} {d["buffer_vol"]:>10.1f}  '
                     f'{d["name"]} ({d["source_conc"]:g} -> {d["final_conc"]:.3g})')
    lines.append('')
    lines.append('  lactaldehyde standard curve (concentrations are IN THE 150 uL reaction)')
    lines.append('  dil.well  from        stock uL  buffer uL  content')
    for g in layout['gradient']:
        origin = (f'dil  {g["source_address"]:<7}' if g['serial']
                  else f'tube {g["source_address"]:<7}')
        lines.append(f'  {g["dilution_well"]:<9} {origin} '
                     f'{g["source_vol"]:>8.1f} {g["buffer_vol"]:>10.1f}  '
                     f'{g["name"]}{"  [serial]" if g["serial"] else ""}')
    return lines


def render_reservoir(layout):
    v = layout['volumes']
    names = {RESERVOIR_BUFFER_WELL: 'assay buffer'}
    names.update({w: 'reaction mix (buffer + NNBT)' for w in RESERVOIR_REACTION_MIX_WELLS})
    names.update({w: 'Purpald reagent' for w in RESERVOIR_PURPALD_WELLS})
    lines = [f'  reservoir ({RESERVOIR_LOADNAME}, slot {SLOT_RESERVOIR})',
             f'  volumes below are what the run consumes - add '
             f'{RESERVOIR_DEAD_VOL_UL / 1000:g} mL dead volume on top of each']
    for well in sorted(layout['reservoir_needs']):
        need = layout['reservoir_needs'][well]
        lines.append(f'    {well}  {names.get(well, "?"):<28} '
                     f'{need / 1000:>5.2f} mL  -> pour '
                     f'{(need + RESERVOIR_DEAD_VOL_UL) / 1000:.1f} mL')
    unused = [w for w in RESERVOIR_REACTION_MIX_WELLS + RESERVOIR_PURPALD_WELLS
              if w not in layout['reservoir_needs']]
    if unused:
        lines.append(f'    {", ".join(unused)}  not needed for this batch - leave empty')
    lines.append('')
    lines.append(f'  tube rack (slot {SLOT_TUBERACK})')
    for d in layout['dilutions']:
        lines.append(f'    {d["source_address"]:<4} {d["name"]:<22} '
                     f'>= {d["source_vol"] + 20:.0f} uL')
    if layout['heat_dilution']:
        lines.append(f'    {TUBE_HEAT_INACT:<4} {HEAT_INACT_NAME:<22} '
                     f'>= {layout["heat_dilution"]["source_vol"] + 20:.0f} uL')
    lines.append(f'    {TUBE_LACTALDEHYDE:<4} {"Lactaldehyde stock":<22} '
                 f'>= {v["lactaldehyde_stock"] + 20:.0f} uL '
                 f'({LAC_STOCK_UM:g} uM)')
    return lines


def deck_contents():
    """{slot: [line, line]} - what goes in each OT-2 deck slot for this protocol."""
    contents = {str(HS_SLOT): ['HEATER-SHAKER', '+ assay plate'],
                SLOT_TUBERACK: ['tube rack', 'enzymes + Lac'],
                SLOT_RESERVOIR: ['reservoir', '12-well 15 mL'],
                SLOT_TIPRACK_20: ['tips 20 uL', ''],
                SLOT_TIPRACK_300: ['tips 300 uL', 'single + 8-ch'],
                SLOT_DILUTION_PLATE: ['dilution plate', '96 flat'],
                '12': ['TRASH', 'fixed, tall']}
    for empty in ('2', '4', '9', '10', '11'):
        contents.setdefault(empty, ['KEEP EMPTY', ''])
    return contents


def render_deck(layout):
    """The OT-2 deck, drawn the way you look at it: slot 1 front-left, 12 back-right."""
    contents = deck_contents()
    w = 18
    bar = '   +' + '+'.join(['-' * w] * 3) + '+'
    lines = ['  OT-2 DECK (slot 1 is front-left, nearest you)', bar]
    for row_start in (10, 7, 4, 1):                     # back row first
        slots = [str(row_start + i) for i in range(3)]
        top, mid, bot = [], [], []
        for s in slots:
            body = contents.get(s, ['', ''])
            top.append(f' {s:<2}'.ljust(w)[:w])
            mid.append(f'  {body[0]}'.ljust(w)[:w])
            bot.append(f'  {body[1]}'.ljust(w)[:w])
        lines += ['   |' + '|'.join(top) + '|',
                  '   |' + '|'.join(mid) + '|',
                  '   |' + '|'.join(bot) + '|', bar]
    lines.append('   Empty slots are REQUIRED, not spare: in single-nozzle mode the')
    lines.append('   pipette overhangs its neighbours. Slot 9 is unusable for')
    lines.append('   single-nozzle work at all - the trash in slot 12 is north of it.')
    lines.append('   See DECK RULES in SECTION 3 before moving anything.')
    return lines


def render_summary(layout):
    v = layout['volumes']
    out = []
    out.append('=' * 78)
    negs = 2 if layout['heat_dilution'] else 1
    out.append(f'  BATCH: {len(ENZYME_BATCH)} enzymes, {len(layout["gradient"])} '
               f'lactaldehyde standards, {negs} negative control(s), all in triplicate '
               f'({REPLICATE_ORIENTATION}-wise)')
    out.append(f'  {layout["wells_used"]} of {PLATE_WELL_COUNT} wells used '
               f'({layout["columns_used"]} columns); room for '
               f'{layout["free_enzyme_slots"]} more enzymes on this plate')
    out.append('=' * 78)
    out.append('')
    out += render_plate_map(layout)
    out.append('')
    out += render_group_key(layout)
    out.append('')
    out.append('-' * 78)
    out += render_dilutions(layout)
    out.append('')
    out.append('-' * 78)
    out += render_reservoir(layout)
    out.append('')
    out.append('-' * 78)
    out += render_deck(layout)
    return out


# =======================================================================================
# SECTION 6 - THE PROTOCOL
# =======================================================================================
def run(protocol):
    # -----------------------------------------------------------------------------
    # Print the whole plan into the run log first. This is the run record: everything
    # you need to read the plate afterwards is in the Opentrons app's log.
    # -----------------------------------------------------------------------------
    for line in render_summary(LAYOUT):
        protocol.comment(line)
    

    # -----------------------------------------------------------------------------
    # DECK + INSTRUMENT SETUP
    # -----------------------------------------------------------------------------
    hs_mod = protocol.load_module('heaterShakerModuleV1', HS_SLOT)
    hs_adapter = hs_mod.load_adapter(HS_ADAPTER_LOADNAME)      # or every Z is undershot
    assay_plate = hs_adapter.load_labware(ASSAY_PLATE_LOADNAME)
    dilution_plate = protocol.load_labware(DILUTION_PLATE_LOADNAME, SLOT_DILUTION_PLATE)
    reservoir = protocol.load_labware(RESERVOIR_LOADNAME, SLOT_RESERVOIR)
    tuberack = protocol.load_labware(TUBERACK_LOADNAME, SLOT_TUBERACK)
    maldi_on = protocol.params.maldi_enabled
    maldi_target = (protocol.load_labware(MALDI_LOADNAME, protocol_api.OFF_DECK)
                    if maldi_on else None)          # swapped onto the deck at U-M-01
    tr_20 = protocol.load_labware('opentrons_96_tiprack_20ul', SLOT_TIPRACK_20)
    tr_300 = protocol.load_labware('opentrons_96_tiprack_300ul', SLOT_TIPRACK_300)

    p20 = protocol.load_instrument('p20_multi_gen2', 'right', tip_racks=[tr_20])
    p300 = protocol.load_instrument('p300_multi_gen2', 'left', tip_racks=[tr_300])

    # pipetting rates
    for pip, rate in ((p20, P20_FLOW_RATE_UL_S), (p300, P300_FLOW_RATE_UL_S)):
        pip.flow_rate.aspirate = rate      # how fast liquid is drawn in
        pip.flow_rate.dispense = rate      # how fast it is pushed out
        p20.flow_rate.blow_out = 7.6           # p20 default; gentle in a shallow well
        p300.flow_rate.blow_out = 94.0         # p300 default - still gentle for this pipette

    # TIPS. Both pipettes now work single-nozzle, and the API refuses starting_tip in a
    # partial configuration, so tips are handed out explicitly. One running index per
    # rack, in rack order (A1, B1, ... H1, A2, ...), shared by single and column pickups
    # so the two can never collide. Set FIRST_TIP_* to start on a part-used rack.
    first_20 = f'{protocol.params.tip20_row}{protocol.params.tip20_col}'    # from the app
    first_300 = f'A{protocol.params.tip300_col}'                            # from the app
    _racks = {'20': tr_20.wells(), '300': tr_300.wells()}
    _next = {'20': _racks['20'].index(tr_20.well(first_20)),
             '300': _racks['300'].index(tr_300.well(first_300))}

    def _take_tip(key, pip, slot, count=1):
        """Hand out `count` consecutive tips; pause for a fresh rack if they run out."""
        if _next[key] + count > len(_racks[key]):
            protocol.pause(f'{key} uL tips used up. Put a FULL rack in slot {slot} '
                           'and resume.')
            _next[key] = 0
        well = _racks[key][_next[key]]                         # row A when count == 8
        _next[key] += count
        pip.pick_up_tip(well)

    def p20_pick_up():
        _take_tip('20', p20, SLOT_TIPRACK_20)                  # one tip

    def p300_pick_up():
        _take_tip('300', p300, SLOT_TIPRACK_300)               # one tip, single-nozzle

    def p300_pick_up_column():
        """8-channel: round up to the next whole column, then take all 8 of it."""
        while _next['300'] % WELLS_PER_COLUMN:                 # skip a part-used column
            _next['300'] += 1
        _take_tip('300', p300, SLOT_TIPRACK_300, count=WELLS_PER_COLUMN)

    protocol.comment(f'tips: 20 uL from {first_20}, 300 uL from {first_300}')


    # ---------------------------------------------------------------------------------
    # ANNOTATION. define_liquid + load_liquid is what makes the Opentrons app draw a
    # coloured, labelled plate map for this run and list the reagent volumes on the
    # "labware setup" screen. It costs nothing at run time and it is the single best
    # way to check the layout is what you meant before you press play.
    # ---------------------------------------------------------------------------------
    liquids = {}
    for region, label in REGION_LABELS.items():
        liquids[region] = protocol.define_liquid(
            name=label, description=f'{label} - {REPLICATES} replicate wells per entry',
            display_color=REGION_COLORS[region])
    for assignment in LAYOUT['assignments']:
        for well_name in assignment['wells']:
            # Annotated with the finished 200 uL well so the map reflects the end state.
            assay_plate[well_name].load_liquid(liquids[assignment['region']], 0)

    rxn_liq = protocol.define_liquid('Reaction mix (buffer + NNBT)', 'premixed off-deck',
                                     '#8e24aa')
    buf_liq = protocol.define_liquid('Assay buffer', 'dilution buffer + blank spike',
                                     '#0288d1')
    pur_liq = protocol.define_liquid('Purpald reagent', 'developing reagent', '#c62828')
    lac_liq = protocol.define_liquid('Lactaldehyde stock',
                                     f'{LAC_STOCK_UM:g} uM positive-control stock',
                                     '#ef6c00')
    enz_liq = protocol.define_liquid('Enzyme stock', 'Bradford-quantified supernatant',
                                     '#2e7d32')

    # Reservoir annotations carry the exact volume to pour, so the app's "labware setup"
    # screen doubles as the prep sheet.
    v = LAYOUT['volumes']
    reservoir_liquids = {RESERVOIR_BUFFER_WELL: buf_liq}
    reservoir_liquids.update({w: rxn_liq for w in RESERVOIR_REACTION_MIX_WELLS})
    reservoir_liquids.update({w: pur_liq for w in RESERVOIR_PURPALD_WELLS})
    for well_name, need in LAYOUT['reservoir_needs'].items():
        reservoir[well_name].load_liquid(reservoir_liquids[well_name],
                                         need + RESERVOIR_DEAD_VOL_UL)
    for d in LAYOUT['dilutions']:
        tuberack[d['source_address']].load_liquid(enz_liq, d['source_vol'] + 20)
    if LAYOUT['heat_dilution']:
        tuberack[TUBE_HEAT_INACT].load_liquid(
            enz_liq, LAYOUT['heat_dilution']['source_vol'] + 20)
    tuberack[TUBE_LACTALDEHYDE].load_liquid(lac_liq, v['lactaldehyde_stock'] + 20)

    # Both pipettes run single-nozzle for the dilution work. Anchored at 'H1' (the
    # front-most nozzle): with 'A1' the gantry cannot reach the assay plate in slot 1.
    p20.configure_nozzle_layout(style=SINGLE, start='H1', tip_racks=[tr_20])
    p300.configure_nozzle_layout(style=SINGLE, start='H1', tip_racks=[tr_300])
    # The p300 is NEVER reconfigured to single-nozzle: that is what lets its rack sit in
    # slot 11 and keeps the middle of the deck free. See DECK RULES in SECTION 3.

    hs_mod.close_labware_latch()

    # ---------------------------------------------------------------------------------
    # SMALL HELPERS
    # ---------------------------------------------------------------------------------
    def pipette_for(volume):
        """p20 up to 20 uL, p300 above it. The p300 doing its own >20 uL transfers is
        the point of this version: a 90 uL move is one trip instead of five."""
        return p20 if volume <= P20_MAX_UL else p300

    def transfer_single(volume, source, dest, pipette=None, mix_after=None,
                        new_tip=True):
        """One p20 transfer, split into 20 uL strokes if it is bigger than the pipette."""
        pip = pipette or pipette_for(volume)
        if new_tip:
            p20_pick_up() if pip is p20 else p300_pick_up()
        remaining = volume
        max_vol = P20_MAX_UL if pip is p20 else P300_MAX_UL
        dispensed = 0.0
        while remaining > 1e-6:
            stroke = min(remaining, max_vol)
            pip.aspirate(stroke, source)
            # Touch off against the SOURCE: any droplet on the outside of the tip falls
            # back into the vessel it caome from rather than riding to the destination.
            touch_off(pip, source)
            pip.dispense(stroke, dest)
            dispensed += stroke
            remaining -= stroke
            if remaining > 1e-6:
                # Mid-transfer: blow out over the liquid already delivered, so every
                # stroke lands its full volume before the next one starts.
                blow_off(pip, dest, above_liquid_ul=dispensed)
        if mix_after:
            reps, mix_vol, well, full_vol = mix_after
            mix_vol = min(mix_vol, max_vol)
            pip.mix(reps, mix_vol, mix_location(well, full_vol, mix_vol))
            pip.blow_out(blowout_location(well, full_vol))
            touch_off(pip, well)
        else:
            # No mix to follow, so this is the last chance to clear the tip.
            blow_off(pip, dest, above_liquid_ul=dispensed)
        if new_tip:
            pip.drop_tip()

    def as_well(target):
        """The Well behind either a Well or a Location (well.bottom(...) etc.)."""
        if hasattr(target, 'labware'):
            return target.labware.as_well()
        return target

    def touch_off(pip, target):
        """Shed the drop on the OUTSIDE of the tip against the well wall."""
        well = as_well(target)
        if TOUCH_TIP_ENABLED and well is not None:
            pip.touch_tip(well, radius=TOUCH_TIP_RADIUS,
                          v_offset=TOUCH_TIP_V_OFFSET_MM)

    def blow_off(pip, target, above_liquid_ul=None):
        """Clear the INSIDE of the tip above any liquid, then touch off.
        `above_liquid_ul` = what the destination already holds; omit for an empty well."""
        well = as_well(target)
        if well is None:
            return
        if BLOW_OUT_ENABLED:
            if above_liquid_ul:
                pip.blow_out(blowout_location(well, above_liquid_ul))
            else:
                pip.blow_out(well.bottom(z=min(BLOWOUT_ABOVE_SURFACE_MM + 1.0,
                                               well.depth - 1.0)))
        touch_off(pip, well)

    def well_area_mm2(well):
        """Cross-sectional area of a well, circular or rectangular."""
        diameter = getattr(well, 'diameter', None)
        if diameter:
            return math.pi * (diameter / 2) ** 2
        return well.length * well.width

    def liquid_height_mm(well, volume_ul):
        return volume_ul / well_area_mm2(well)

    def mix_location(well, full_vol_ul, stroke_vol_ul):
        """Tip height that stays submerged for the WHOLE mix stroke. The level that
        matters is the one at the BOTTOM of the stroke, not at rest."""
        low = liquid_height_mm(well, max(0.0, full_vol_ul - stroke_vol_ul))
        z = max(DILUTION_MIX_MIN_HEIGHT_MM, low * DILUTION_MIX_DEPTH_FRACTION)
        return well.bottom(z=min(z, well.depth - 1.0))

    def blowout_location(well, full_vol_ul):
        """Just above the surface: blowing out under it foams a protein solution."""
        z = liquid_height_mm(well, full_vol_ul) + BLOWOUT_ABOVE_SURFACE_MM
        return well.bottom(z=min(z, well.depth - 1.0))

    def aspirate_depth(well, volume_left_ul, stroke_ul):
        """Tip height that stays submerged for the WHOLE aspirate stroke - the level
        AFTER it, not before. Every fixed aspirate height here has been wrong once:
        2.0 mm missed the gradient source, the API's 1.0 mm default ran dry on the
        third replicate of a spike. So it is computed wherever the volume is known."""
        after = max(0.0, volume_left_ul - stroke_ul)
        z = max(DILUTION_ASPIRATE_MIN_HEIGHT_MM, liquid_height_mm(well, after) * DILUTION_ASPIRATE_DEPTH_FRACTION)
        return well.bottom(z=min(z, well.depth - 1.0))

    def dil(well_name):
        return dilution_plate[well_name].bottom(z=DILUTION_DISPENSE_HEIGHT_MM)

    def load_triplicate(source_location, well_names, source_volume_ul=None):
        """One spike into each replicate well, ONE FRESH ASPIRATE EACH (not one big
        aspirate split three ways - that is what keeps replicates comparable).
        One tip per group. `source_volume_ul` = what the source holds, so the tip can
        follow the surface down; None for a reservoir, which never runs low."""
        p20_pick_up()
        left = source_volume_ul
        for name in well_names:
            # Follow the surface down: three 10 uL draws take a 60 uL well to 30 uL,
            # and a fixed 1 mm tip is above the liquid by the third one.
            src = (aspirate_depth(as_well(source_location), left, SPIKE_VOL_UL)
                   if source_volume_ul else source_location)
            if left is not None:
                left -= SPIKE_VOL_UL
            p20.aspirate(SPIKE_VOL_UL, src)
            # Touch off against the source so nothing rides across on the outside.
            touch_off(p20, source_location)
            p20.dispense(SPIKE_VOL_UL,
                         assay_plate[name].bottom(z=SPIKE_DISPENSE_HEIGHT_MM))
            # 10 uL released 6 mm above a DRY well does not always fall on its own -
            # surface tension holds it on the tip. Blow out and touch off so the full
            # spike actually lands, and every replicate gets the same volume.
            blow_off(p20, assay_plate[name])
        p20.drop_tip()

    all_dilutions = LAYOUT['all_dilutions']
    mix_vol = DILUTION_WELL_VOL_UL * DILUTION_MIX_STROKE_FRACTION

    # --- U-O-00a  buffer into every dilution well --------------------------------
    # Empty wells, so one tip serves the whole step - nothing it touches can be
    # contaminated. The tip is picked up lazily, the first time it is needed.
    protocol.comment('U-O-00a  filling dilution wells with buffer')
    tipped = set()
    for entry in all_dilutions + LAYOUT['gradient']:
        if entry['buffer_vol'] <= 0:
            continue
        pip = pipette_for(entry['buffer_vol'])
        if pip not in tipped:
            p20_pick_up() if pip is p20 else p300_pick_up()
            tipped.add(pip)
        transfer_single(entry['buffer_vol'], reservoir[RESERVOIR_BUFFER_WELL],
                        dil(entry['dilution_well']), pipette=pip, new_tip=False)
    for pip in tipped:
        pip.drop_tip()

    # --- U-O-00b  enzyme stock -> dilution well, then mix -------------------------
    # Fresh tip per enzyme: different physical samples, and the tip goes INTO the well
    # to mix, so a reused one would carry enzyme across.
    _probe = dilution_plate['A1']
    _stroke = min(mix_vol, P20_MAX_UL)
    protocol.comment(
        f'    mixing: {DILUTION_WELL_VOL_UL:g} uL stands '
        f'{liquid_height_mm(_probe, DILUTION_WELL_VOL_UL):.2f} mm deep; '
        f'{_stroke:g} uL strokes bottom out at '
        f'{liquid_height_mm(_probe, DILUTION_WELL_VOL_UL - _stroke):.2f} mm; '
        f'tip sits at '
        f'{mix_location(_probe, DILUTION_WELL_VOL_UL, _stroke).point.z - _probe.bottom().point.z:.2f} mm')
    protocol.comment('U-O-00b  diluting enzymes to a common concentration')
    for d in all_dilutions:
        protocol.comment(
            f'    {d["name"]}: {d["source_vol"]:.1f} uL from tube {d["source_address"]} '
            f'+ {d["buffer_vol"]:.1f} uL buffer -> {d["dilution_well"]}')
        transfer_single(d['source_vol'], tuberack[d['source_address']],
                        dil(d['dilution_well']),
                        mix_after=(DILUTION_MIX_REPS, mix_vol,
                                   dilution_plate[d['dilution_well']],
                                   DILUTION_WELL_VOL_UL))

    # --- U-O-00c  lactaldehyde standard curve ------------------------------------
    # Most concentrated point first, so a [serial] source is always finished before
    # anything draws from it. Fresh tip per point.
    protocol.comment('U-O-00c  building the lactaldehyde standard curve')
    for g in LAYOUT['gradient']:
        src = (aspirate_depth(dilution_plate[g['source_address']],
                              DILUTION_WELL_VOL_UL, g['source_vol'])
               if g['serial'] else tuberack[g['source_address']])
        protocol.comment(
            f'    {g["name"]}: {g["source_vol"]:.1f} uL from '
            f'{"dilution well" if g["serial"] else "tube"} {g["source_address"]} '
            f'+ {g["buffer_vol"]:.1f} uL buffer -> {g["dilution_well"]}')
        transfer_single(g['source_vol'], src, dil(g['dilution_well']),
                        mix_after=(DILUTION_MIX_REPS, mix_vol,
                                   dilution_plate[g['dilution_well']],
                                   DILUTION_WELL_VOL_UL))

    # --- U-O-01  spikes into the assay plate -------------------------------------
    # One loop covers enzymes, both negative controls and the whole curve: every
    # source was already resolved in build_layout(), in plate-map order.
    protocol.comment('U-O-01..04  loading 10 uL spikes into the assay plate')
    for a in LAYOUT['assignments']:
        kind, address = a['source']
        # A reservoir well is deep and effectively bottomless for a 10 uL draw; a
        # dilution well holds DILUTION_WELL_VOL_UL and visibly drains as we take from it.
        source = (reservoir[address] if kind == 'reservoir'
                  else dilution_plate[address])
        source_volume = None if kind == 'reservoir' else DILUTION_WELL_VOL_UL
        protocol.comment(f'    {a["label"]} -> {", ".join(a["wells"])}')
        load_triplicate(source, a['wells'], source_volume)

    # --- U-O-02  reaction mix, 8-channel, whole columns --------------------------
    # One aspirate per dispense: a shared aspirate leaves an air pocket between
    # dispenses, which is what lets droplets fall between columns. Flow rate halved
    # for this step only - the NNBT is in acetonitrile and drips at full speed.
    protocol.comment('U-O-05  adding reaction mix')
    p300.configure_nozzle_layout(style=ALL, tip_racks=[tr_300])   # back to 8-channel
    plate_columns = [assay_plate[w].bottom(z=REAGENT_DISPENSE_HEIGHT_MM)
                     for w in LAYOUT['plate_columns']]

    default_aspirate, default_dispense = p300.flow_rate.aspirate, p300.flow_rate.dispense
    p300.flow_rate.aspirate = default_aspirate * REACTION_MIX_FLOW_RATE_SCALE
    p300.flow_rate.dispense = default_dispense * REACTION_MIX_FLOW_RATE_SCALE
    p300_pick_up_column()
    for column, src_well in zip(plate_columns, LAYOUT['rxn_well_per_column']):
        p300.aspirate(REACTION_MIX_VOL_UL, reservoir[src_well])
        # Touch off in the reservoir: the mix is the viscous one, and a droplet on the
        # outside of eight tips is 8 x a droplet of loss per column.
        if TOUCH_TIP_ENABLED:
            p300.touch_tip(reservoir[src_well], radius=TOUCH_TIP_RADIUS,
                           v_offset=TOUCH_TIP_V_OFFSET_MM)
        p300.dispense(REACTION_MIX_VOL_UL, column)
        # Blow out at the dispense height, which is well above the 150 uL surface, so
        # the full 140 uL lands and nothing rides to the next column.
        if BLOW_OUT_ENABLED:
            p300.blow_out(column)
        p300.touch_tip()
    p300.drop_tip()
    p300.flow_rate.aspirate = default_aspirate
    p300.flow_rate.dispense = default_dispense

    # --- U-M-01  PAUSE - parafilm on (2 h at 40 C evaporates an open plate) ------
    # With MALDI on the plate CANNOT be covered: the robot reaches into it every round.
    # It will evaporate, which concentrates what is left and biases the course upward.
    if maldi_on:
        protocol.move_labware(tuberack, protocol_api.OFF_DECK, use_gripper=False)
        protocol.move_labware(maldi_target, MALDI_SLOT, use_gripper=False)
        protocol.pause(
            f'MALDI time course: leave the plate UNCOVERED and seat the target in slot '
            f'{MALDI_SLOT}. It will evaporate over {INCUBATION_MIN} min - expected.')
    else:
        protocol.pause(
            f'Cover the assay plate with parafilm, then resume for the {INCUBATION_MIN} '
            f'min incubation at {INCUBATION_TEMP_C} degC.')

    # --- U-O-03  incubate --------------------------------------------------------
    # Heater stays ON afterwards: development runs at the same temperature.
    if INCUBATION_TEMP_C >= 37:
        hs_mod.set_and_wait_for_temperature(INCUBATION_TEMP_C)
    else:
        protocol.comment(
            f'INCUBATION_TEMP_C ({INCUBATION_TEMP_C}) is below the heater-shaker\'s '
            '37-95 degC controllable range: shaking at ambient instead.')

    if not maldi_on:
        hs_mod.set_and_wait_for_shake_speed(INCUBATION_RPM)
        protocol.delay(minutes=INCUBATION_MIN)
        hs_mod.deactivate_shaker()
    else:
        # Spot at t = 0, then every interval. Each round gets a fresh block of rows, so
        # the target reads as one row-block per timepoint, in plate-map order.
        interval = protocol.params.maldi_interval
        rounds = int(INCUBATION_MIN // interval) + 1
        per_row = MALDI_COLUMNS // MALDI_SPOT_REPLICATES
        rows_per_round = math.ceil(len(LAYOUT['assignments']) / per_row)
        if rows_per_round * rounds > len(MALDI_ROWS):
            raise ValueError(
                f'{len(LAYOUT["assignments"])} conditions x {rounds} rounds needs '
                f'{rows_per_round * rounds} target rows, only {len(MALDI_ROWS)} exist. '
                'Raise the interval or spot fewer conditions.')
        draw = MALDI_SPOT_VOL_UL * MALDI_SPOT_REPLICATES + MALDI_ASPIRATE_EXTRA_UL

        for r in range(rounds):
            if r:                                   # shake between rounds, not during
                hs_mod.set_and_wait_for_shake_speed(INCUBATION_RPM)
                protocol.delay(minutes=interval)
                hs_mod.deactivate_shaker()          # never pipette a moving plate
            protocol.comment(f'MALDI round {r + 1}/{rounds}, nominal t = '
                             f'{r * interval} min')
            for i, a in enumerate(LAYOUT['assignments']):
                row = MALDI_ROWS[r * rows_per_round + i // per_row]
                col0 = (i % per_row) * MALDI_SPOT_REPLICATES + 1
                spots = [f'{row}{col0 + k}' for k in range(MALDI_SPOT_REPLICATES)]
                p20_pick_up()                       # fresh tip: plate then target
                p20.aspirate(draw, assay_plate[a['wells'][0]].bottom(
                    z=MALDI_SAMPLE_HEIGHT_MM))      # always well 1; 2 and 3 stay clean
                for sp in spots:
                    p20.dispense(MALDI_SPOT_VOL_UL,
                                 maldi_target[sp].bottom(z=MALDI_SPOT_HEIGHT_MM))
                p20.drop_tip()                      # surplus leaves with the tip
                protocol.comment(f'    {a["label"]}: {a["wells"][0]} -> '
                                 f'{", ".join(spots)}')
        leftover = INCUBATION_MIN - (rounds - 1) * interval
        if leftover > 0:                            # finish the incubation
            hs_mod.set_and_wait_for_shake_speed(INCUBATION_RPM)
            protocol.delay(minutes=leftover)
            hs_mod.deactivate_shaker()
        protocol.pause(f'Take the MALDI target from slot {MALDI_SLOT}, dry it and apply '
                       'matrix. Spot key is in this run log. Resume to finish.')

    # --- U-M-02  PAUSE - parafilm off --------------------------------------------
    protocol.pause('Remove the parafilm from the assay plate, then resume.')

    # --- U-O-04  Purpald -> 200 uL -----------------------------------------------
    # Shaker already stopped: dispensing into a moving plate misplaces the droplet.
    protocol.comment('U-O-07  adding Purpald')
    p300_pick_up_column()
    for column, src_well in zip(plate_columns, LAYOUT['purpald_well_per_column']):
        p300.aspirate(PURPALD_VOL_UL, reservoir[src_well])
        if TOUCH_TIP_ENABLED:
            p300.touch_tip(reservoir[src_well], radius=TOUCH_TIP_RADIUS,
                           v_offset=TOUCH_TIP_V_OFFSET_MM)
        p300.dispense(PURPALD_VOL_UL, column)
        if BLOW_OUT_ENABLED:
            p300.blow_out(column)
        p300.touch_tip()
    p300.drop_tip()

    # --- U-O-05  develop ---------------------------------------------------------
    # 1000 rpm, not the module's 3000: parafilm off and wells brim-full, and a
    # slopping plate cross-contaminates a triplicate assay.
    if DEVELOP_TEMP_C >= 37:
        hs_mod.set_and_wait_for_temperature(DEVELOP_TEMP_C)
    hs_mod.set_and_wait_for_shake_speed(DEVELOP_RPM)
    protocol.delay(minutes=DEVELOP_MIN)
    hs_mod.deactivate_shaker()
    hs_mod.deactivate_heater()

    # --- U-C-06  hand off to the plate reader (A530) -----------------------------
    hs_mod.open_labware_latch()
    protocol.comment('Done - read A530. Group key is printed at the top of this log.')


# =======================================================================================
# SECTION 7 - LOCAL PREVIEW
# `python3 <this file>` needs no robot and no opentrons install: prints the plate map,
# group key, dilution recipes, deck map and what to load, and writes plate_map.html.
# =======================================================================================
HTML_TEMPLATE = """<!doctype html>
<meta charset="utf-8"><title>{title}</title>
<style>
 body{{font:14px/1.5 system-ui,-apple-system,Segoe UI,sans-serif;margin:24px;color:#1a1a1a}}
 h1{{font-size:20px;margin:0 0 4px}} h2{{font-size:15px;margin:28px 0 8px}}
 .sub{{color:#666;margin-bottom:20px}}
 table{{border-collapse:separate;border-spacing:3px}}
 th{{font-weight:600;color:#666;font-size:12px}}
 td.w{{width:66px;height:52px;border-radius:6px;text-align:center;vertical-align:middle;
   color:#fff;font-size:11px;line-height:1.25;padding:2px;overflow:hidden}}
 td.empty{{background:#f0f0f0;color:#bbb}}
 .n{{font-weight:700;font-size:12px;display:block}}
 ul.legend{{list-style:none;padding:0}} ul.legend li{{margin:4px 0}}
 .sw{{display:inline-block;width:14px;height:14px;border-radius:3px;
   vertical-align:-2px;margin-right:8px}}
 pre{{background:#f7f7f7;padding:14px;border-radius:8px;overflow:auto;font-size:12px}}
</style>
<h1>{title}</h1>
<div class="sub">{subtitle}</div>
<table>{table}</table>
<h2>Legend</h2><ul class="legend">{legend}</ul>
<h2>Group key</h2><pre>{groups}</pre>
<h2>Dilution recipes</h2><pre>{dilutions}</pre>
<h2>What to load</h2><pre>{reservoir}</pre>
<h2>Deck</h2><pre>{deck}</pre>
"""


def write_html(layout, path='plate_map.html'):
    index_of = {}
    for i, a in enumerate(layout['assignments']):
        for w in a['wells']:
            index_of[w] = i + 1
    header = '<tr><th></th>' + ''.join(f'<th>{c}</th>' for c in range(1, 13)) + '</tr>'
    rows = [header]
    for row in ROW_LETTERS:
        cells = [f'<th>{row}</th>']
        for col in range(1, 13):
            w = f'{row}{col}'
            info = layout['plate_map'].get(w)
            if not info:
                cells.append(f'<td class="w empty">{w}</td>')
            else:
                cells.append(
                    f'<td class="w" style="background:{REGION_COLORS[info["region"]]}">'
                    f'<span class="n">{w}</span>{info["label"]}<br>rep {info["replicate"]}'
                    '</td>')
        rows.append('<tr>' + ''.join(cells) + '</tr>')
    legend = ''.join(
        f'<li><span class="sw" style="background:{REGION_COLORS[r]}"></span>{l}</li>'
        for r, l in REGION_LABELS.items())
    html = HTML_TEMPLATE.format(
        title='W-O-SC-03 plate map',
        subtitle=(f'{len(ENZYME_BATCH)} enzymes &middot; {len(layout["gradient"])} '
                  f'lactaldehyde standards &middot; '
                  f'{2 if layout["heat_dilution"] else 1} negative control(s) &middot; '
                  f'{REPLICATE_ORIENTATION}-wise triplicates &middot; '
                  f'{layout["wells_used"]}/96 wells &middot; '
                  f'room for {layout["free_enzyme_slots"]} more enzymes'),
        table=''.join(rows), legend=legend,
        groups='\n'.join(render_group_key(layout)),
        dilutions='\n'.join(render_dilutions(layout)),
        reservoir='\n'.join(render_reservoir(layout)),
        deck='\n'.join(render_deck(layout)))
    with open(path, 'w') as fh:
        fh.write(html)
    return path


if __name__ == '__main__':
    print('\n'.join(render_summary(LAYOUT)))
    print()
    print(f'wrote {write_html(LAYOUT)}')
