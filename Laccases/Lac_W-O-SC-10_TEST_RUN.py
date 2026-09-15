"""
W-O-SC-10   NNBT/Purpald + guaiacol + ABTS + MALDI       Kristine Toft Johansen s215098
=======================================================================================
One dilution series feeds three readouts. Everything you change per run is a RUNTIME
PARAMETER in the Opentrons app - you should not need to open this file.

  NNBT plate  (slot 1, Heater-Shaker)  two arms: without and with guaiacol
  ABTS plate  (slot 5, swapped in)     activity, read immediately at A414/A734
  MALDI target(slot 5, swapped again)  spotted during the NNBT incubation

NNBT PLATE - 32 triplicate groups, 8 per 3-column block
  cols  1-3   NC1, NC2, NC3, 5 lactaldehyde standards        mix: NNBT (NC3: no-NNBT)
  cols  4-6   the same 8, with guaiacol                      mix: guaiacol (NC3: no-NNBT+gua)
  cols  7-9   enzymes 1-8                                    mix: NNBT
  cols 10-12  enzymes 1-8, with guaiacol                     mix: guaiacol
  -> 8 enzymes maximum. NC1 = buffer only, NC2 = heat-inactivated enzyme,
     NC3 = active enzyme in reaction mix WITHOUT NNBT.

WELL RECIPE   10 uL spike + 140 uL reaction mix = 150 uL, 2 h @ 40 C,
              + 50 uL Purpald = 200 uL (the plate's full capacity)

MULTI-DISPENSE   wherever one reagent goes to several wells in a row, ONE aspirate
serves several of them, with a touch-tip after every dispense: spikes 2 wells per
p20 load, reaction mix 2 wells per p300 load, Purpald 5 columns per load. Each p300
load carries a 20 uL disposal volume that is blown into the trash - build_layout()
has already added that to the volumes it asks you to pour. TWO PLACES DO NOT DO THIS:
the ABTS mix, because 2 x 190 uL will not fit a p300, and everything on the MALDI
target, where a touch-tip is forbidden (see the spotting section).

SLOT 5 IS USED THREE TIMES   tube rack -> ABTS plate -> MALDI target. They are never
needed at once. It is the middle deck column, so BOTH pipettes can reach every well of
everything that lands there - see rule 3 below.

RUN ORDER
  1  buffer into the dilution wells
  2  enzyme / NC3 / heat-inactivated stocks in, mixed      -> all at one molar conc
  3  lactaldehyde standard curve built by serial dilution
  4  ABTS standards moved from tubes into dilution wells   (rack is about to leave)
  5  10 uL spikes into the NNBT plate, all four blocks
  6  PAUSE - swap the tube rack for the ABTS plate
  7  10 uL spikes into the ABTS plate
  8  reaction mixes into the NNBT plate
       blocks 1-2 one well at a time  (NC3 needs a different mix from its neighbours)
       blocks 3-4 8-channel by column (all one mix)
  9  ABTS mix - LAST, because ABTS is kinetic and starts on contact
 10  PAUSE - ABTS plate to the reader, seal the NNBT plate, MALDI labware in
 11  incubate 2 h @ 40 C, pausing every interval to unseal, spot MALDI, reseal
 12  Purpald, develop, read A530

TWO RULES THAT HAVE ALREADY COST A RUN
  1  Slot 9 can never take a single-nozzle move: the trash in slot 12 is north of it,
     and opentrons_simulate does not model that. check_deck() enforces it.
  2  A HEIGHT is not a DEPTH. A tip above the liquid draws air; a tip on a flat well
     bottom seals and moves nothing. Anything that must reach liquid is COMPUTED.
  3  The LEFT mount cannot reach the right-hand deck column. Both mounts ride one X
     carriage ~34 mm apart, against an X hard limit at 418 mm, so the left pipette
     (p300) tops out near deck X 382. Slots 3/6/9 start at x=265, so a well past
     x~117 in one of them is unreachable - and touch_tip adds another 4 mm. That is
     what killed the 2026-09-11 run: the tube rack was in slot 3 and a touch_tip on
     tube D6 commanded X422.2 -> 'ALARM: Hard limit +X'. SLOT 3 IS NOW FOR THE 20 uL
     TIP RACK ONLY, which only the right-mount p20 ever visits. Nothing the p300
     touches may go in the right-hand column unless its wells stay left of x~113.
"""

import math

try:                                          # only present on the robot / simulator
    from opentrons import protocol_api
    from opentrons.protocol_api import SINGLE, ALL
except ImportError:                           # pragma: no cover - local preview
    protocol_api = None
    SINGLE = ALL = None

MAX_ENZYMES = 8                               # 4 blocks of 8 groups fills the plate


def add_parameters(p):
    """Everything you set per run. Shown in the app at run setup."""
    p.add_int(display_name='Number of enzymes', variable_name='n_enz', default=1,
              minimum=1, maximum=MAX_ENZYMES,
              description=f'TEST BUILD default 1. Max {MAX_ENZYMES}.')
    p.add_float(display_name='Target conc (uM)', variable_name='target_um',
                default=0.0, minimum=0.0, maximum=1000.0,
                description='Molar conc every enzyme is diluted to. 0 = auto.')
    for i in range(1, MAX_ENZYMES + 1):       # only the first n_enz are used
        p.add_float(display_name=f'{i}. Bradford (mg/mL)', variable_name=f'mg_ml_{i}',
                    default=1.00, minimum=0.0, maximum=100.0,
                    description=f'Enzyme {i}, tube {_tube(i - 1)}. Ignored above n.')
        p.add_float(display_name=f'{i}. Molar mass (kDa)', variable_name=f'mw_{i}',
                    default=65.0, minimum=1.0, maximum=500.0,
                    description=f'Enzyme {i} molar mass. Ignored above n.')
    p.add_float(display_name='NC2 heat-inact (mg/mL)', variable_name='heat_mg_ml',
                default=2.40, minimum=0.01, maximum=100.0,
                description='Heat-inactivated control, tube C6.')
    p.add_float(display_name='NC2 heat-inact (kDa)', variable_name='heat_mw',
                default=65.0, minimum=1.0, maximum=500.0,
                description='Heat-inactivated control molar mass.')
    p.add_int(display_name='NC3 uses enzyme #', variable_name='nc3_enz', default=1,
              minimum=1, maximum=MAX_ENZYMES,
              description='Which enzyme goes in NC3 (mix without NNBT).')
    # --- tips ------------------------------------------------------------------------
    p.add_str(display_name='20 uL tip row', variable_name='tip20_row', default='A',
              choices=[{'display_name': r, 'value': r} for r in 'ABCDEFGH'],
              description='Row of the first UNUSED tip in the 20 uL rack.')
    p.add_int(display_name='20 uL tip column', variable_name='tip20_col', default=1,
              minimum=1, maximum=12, description='Column of that first unused tip.')
    p.add_int(display_name='300 uL tip column', variable_name='tip300_col', default=1,
              minimum=1, maximum=12,
              description='First UNUSED column in the 300 uL rack.')
    # --- MALDI -----------------------------------------------------------------------
    p.add_bool(display_name='MALDI spotting', variable_name='maldi_on', default=True,
               description='TEST BUILD default ON, so spotting gets tested.')
    p.add_int(display_name='MALDI interval (min)', variable_name='maldi_interval',
              default=30, minimum=5, maximum=120,
              description='Minutes between spotting rounds.')
    p.add_str(display_name='MALDI start row', variable_name='maldi_row', default='A',
              choices=[{'display_name': r, 'value': r} for r in 'ABCDEFGHIJKLMNOP'],
              description='First target row. Each round uses two rows.')


metadata = {
    'apiLevel': '2.20',
    'protocolName': 'W-O-SC-10 TEST RUN - one pass through every step',
    'description': (
        'SHORTENED TEST BUILD - NOT THE ASSAY. One enzyme, duplicates instead of '
        'triplicates, 2 min incubate and 2 min develop instead of 120 and 10. Both '
        'standard curves are full length. Every mechanism the real run uses is '
        'exercised at least once: both pipettes, mixing, multi-dispense pairing, '
        'touch-tip and blow-out, the Heater-Shaker (latch, heat, shake), all three '
        'off-deck plate swaps, and one round of MALDI spotting. Search TEST BUILD in '
        'the file for every line that differs from the assay.'),
}


# =======================================================================================
# SECTION 1 - THINGS THAT RARELY CHANGE   (everything per-run is a parameter above)
# =======================================================================================
ENZYME_NAMES = ['Lac-01', 'Lac-02', 'Lac-03', 'Lac-04',
                'Lac-05', 'Lac-06', 'Lac-07', 'Lac-08']     # plate-map labels
ENZ_TARGET_AUTO_FRACTION = 0.90        # auto target = this x the weakest enzyme, in uM

# Lactaldehyde positive control. FIVE points - the sixth was dropped to make room for
# NC3, so the controls fill one 3-column block exactly.
LAC_STOCK_UM = 100_000.0               # tube D6. 1 M diluted 1:10 = 100000 uM
LAC_GRADIENT_UM = [500, 1500, 3000, 4500, 6000]   # uM IN THE 150 uL REACTION
LAC_SERIAL = False                     # False: every standard comes straight from the
                                       # D6 stock. True: the old chain, each point made
                                       # from the one above it.
                                       # The chain CANNOT WORK with this gradient. The
                                       # 6000 point has to hand 3/4 of itself to the
                                       # 4500 point and then still give 60 uL to its own
                                       # spikes - that is 210 uL out of a 200 uL well.
                                       # It ran dry, and the p300 was drawing air off
                                       # the bottom of an empty well. Direct dilution
                                       # also stops dilution error compounding down the
                                       # chain. Set True only if you raise
                                       # DILUTION_WELL_VOL_UL past ~400 uL, which this
                                       # plate cannot hold.

# ABTS protein standards - PRE-MADE BY HAND, one tube each, in this order.
ABTS_STD_MG_ML = [0.0, 0.125, 0.25, 0.5, 0.75, 1.0]
ABTS_STD_TUBES = ['A5', 'B5', 'C5', 'D5', 'A6', 'B6']

DILUTION_WELL_VOL_UL = 200.0           # per dilution well. An enzyme well gives up 90 uL
                                       # (NNBT arm A + arm B + ABTS), and a serial lac
                                       # step gives up as much as 3/4 of the well. At
                                       # 150 that left the F1<-E1 draw with 0.22 mm of
                                       # liquid over the tip - the 0.8 mm minimum height
                                       # overrode the half-depth rule and the p300 drew
                                       # air. The RATIOS are unchanged, so no
                                       # concentration changes; there is simply more
                                       # left behind. 200 restores 0.56 mm on that draw.


# =======================================================================================
# SECTION 2 - ASSAY CHEMISTRY
# =======================================================================================
SPIKE_VOL_UL = 10.0                    # sample per NNBT / ABTS well
REACTION_MIX_VOL_UL = 140.0            # 10 + 140 = 150 uL during the incubation
PURPALD_VOL_UL = 50.0                  # -> 200 uL, the plate's capacity
ABTS_MIX_VOL_UL = 190.0                # 10 + 190 = 200 uL on the ABTS plate
REACTION_VOL_UL = SPIKE_VOL_UL + REACTION_MIX_VOL_UL

REPLICATES = 2                         # TEST BUILD: 2, not 3. Not 1 - multi_dispense
                                       # pairs TWO destinations onto one p20 aspirate
                                       # and gives the blow-out to the second, so a
                                       # single replicate would never run that path.
                                       # 2 is the smallest number that still tests it.
WELLS_PER_COLUMN = 8
ROW_LETTERS = list('ABCDEFGH')
BLOCK_COLS = REPLICATES                # a block is REPLICATES columns wide

INCUBATION_TEMP_C = 40                 # inside the module's 37-95 C range
INCUBATION_RPM = 250
INCUBATION_MIN = 2                     # TEST BUILD: 2, assay value is 120. Also drops
                                       # MALDI to a single spotting round, since
                                       # rounds = INCUBATION_MIN // interval + 1.
DEVELOP_TEMP_C = 40
DEVELOP_RPM = 1000                     # not 3000: parafilm off, wells brim-full
DEVELOP_MIN = 2                        # TEST BUILD: 2, assay value is 10
REACTION_MIX_FLOW_SCALE = 0.5          # NNBT is in acetonitrile and drips at full speed
# Pipetting speed. The defaults (7.6 / 94 uL/s) make mixing crawl; these are your
# SC-08 values. Blow-out stays at each pipette's own default - a fast blow-out into a
# shallow well spatters and foams protein.
P20_FLOW_UL_S = 20.0                   # default 7.6, hardware max 24
P300_FLOW_UL_S = 200.0                 # default 94, hardware max 275
P20_BLOWOUT_UL_S = 10                 # p20 default: gentle
P300_BLOWOUT_UL_S = 110.0               # p300 default: gentle for this pipette
MALDI_FLOW_UL_S = 7.6                  # spotting only. A 1 uL droplet placed onto a
                                       # flat steel target at 15 uL/s splashes.

# MALDI: 5 uL sample + 20 uL milliQ = 1:5, made in a fresh plate, one spot per condition.
MALDI_SAMPLE_UL = 5.0
MALDI_WATER_UL = 20.0
MALDI_SPOT_UL = 1.0                    # sample per spot; p20 minimum
MALDI_MATRIX_UL = 1.0                  # matrix goes down FIRST, sample lands on it
MALDI_MIX_REPS = 3                     # mixed in place on the target (Lukas W-O-MS-01)
MALDI_MIX_UL = 1.0                     # the p20's minimum, and half the 2 uL spot
# MATRIX DRIES. Matrix is laid a batch at a time, then that batch's samples go on
# immediately - so no matrix drop waits longer than one batch. Smaller batch = shorter
# wait but more tips (one extra matrix tip per batch); bigger = fewer tips, longer wait.
MALDI_MATRIX_BATCH = 4
MALDI_SPOT_HEIGHT_MM = 0.3             # from the WELL BOTTOM, and the well is only
                                       # 0.1 mm deep, so this is ~0.2 mm above the
                                       # target face: into the droplet, not pressed on
                                       # the steel (which would occlude the orifice).
                                       # The API's ~1 mm default would sit 0.9 mm clear
                                       # and the drop would FALL instead of be placed.
                                       # Lukas's W-O-MS-01 value, matched to this
                                       # labware definition. DRY-RUN IT ON A BARE TARGET.
MALDI_DRAW_HEIGHT_MM = 3.0             # draw height in the 150 uL NNBT well
MALDI_ROWS = list('ABCDEFGHIJKLMNOP')  # 16 rows on the target
MALDI_COLS = 24


# =======================================================================================
# SECTION 3 - HARDWARE
# =======================================================================================
# DECK. In single-nozzle mode the pipette body overhangs the slot to the NORTH (slot+3),
# so nothing tall may sit there. Tall = tip rack (64 mm), tube rack (54 mm); a plate
# (14 mm) or the reservoir (31 mm) is fine. check_deck() enforces it.
HS_SLOT = 1                            # Heater-Shaker + NNBT plate  (north: 4, empty)
SLOT_SWAP = '5'                        # tube rack -> ABTS plate -> MALDI target
                                       #                             (north: 8, reservoir)
SLOT_TIPRACK_20 = '3'                  # 20 uL tips, p20 ONLY - see rule 3 in the header
                                       #                             (north: 6, dil plate)
SLOT_DILUTION = '6'                    # dilution plate -> MALDI dilution plate
SLOT_TIPRACK_300 = '7'                 # 300 uL tips                 (north: 10, empty)
SLOT_RESERVOIR = '8'                   # reagents                    (north: 11, empty)
# Slots 2, 4, 9, 10, 11 MUST STAY EMPTY - clearance, not spare space.

NNBT_PLATE = 'eppendorf_96_wellplate_350ul'
DILUTION_PLATE = 'corning_96_wellplate_360ul_flat'
ABTS_PLATE = 'corning_96_wellplate_360ul_flat'
MALDI_PLATE = 'maldi_384_wellplate'    # Lukas's real definition, face at 18.0 mm
RESERVOIR_LOADNAME = 'nest_12_reservoir_15ml'
TUBERACK = '3d_printed_tuberack_1.5ml'
HS_ADAPTER = 'opentrons_96_flat_bottom_adapter'  # MUST be declared or every Z is wrong

TUBERACK_ROWS = 4                      # rack is 4 rows (A-D) x 6 columns = 24 tubes
TUBE_HEAT_INACT = 'C6'                 # fixed positions, so they never move
TUBE_LACTALDEHYDE = 'D6'

# Reservoir. Four different reaction mixes, all premixed off-deck by you.
RES_NNBT = ['A1', 'A2']                # buffer + NNBT                (2 wells: >13 mL)
RES_BUFFER = 'A3'                      # dilution buffer + the NC1 spike
RES_PURPALD = ['A4']
RES_ABTS = ['A5', 'A6']                # ABTS reaction mix
RES_GUAIACOL = ['A7', 'A8']            # guaiacol buffer + NNBT + guaiacol
RES_NO_NNBT = ['A9']                   # NC3: reaction mix WITHOUT NNBT
RES_NO_NNBT_GUA = ['A10']              # NC3 guaiacol: no NNBT, with guaiacol
RES_WATER = 'A11'                      # milliQ for the MALDI 1:5 dilution
RES_MATRIX = 'A12'                     # MALDI matrix - the robot spots it
RES_USABLE_UL = 13_000.0               # 15 mL nominal, minus fill margin
RES_DEAD_UL = 1_000.0                  # pour this much extra so tips never hit air

P20_MIN, P20_MAX = 1.0, 20.0
P300_MIN, P300_MAX = 20.0, 300.0

# --- tip heights: absolute, above the well bottom, for staying OUT of liquid ----------
EXTRA_HEIGHT_MM = 2.0                  # padding for the residual real-vs-modelled Z gap
REAGENT_HEIGHT_MM = 13.0 + EXTRA_HEIGHT_MM   # reagents dropped onto liquid
SPIKE_HEIGHT_MM = 1.0 + EXTRA_HEIGHT_MM      # spikes into empty wells
SPIKE_LOW_HEIGHT_MM = 1.5              # a multi-dispensed spike gets NO blow-out (the
                                       # tip still holds the next well's 10 uL), so it
                                       # is placed low enough for the droplet to meet
                                       # the floor and be shed by the touch-tip
DILUTION_HEIGHT_MM = 0.5               # dispensing into a dilution well
# --- computed depths: for reaching INTO liquid ---------------------------------------
MIX_REPS = 3
MIX_STROKE_FRACTION = 0.4              # of the well volume
MIX_MIN_HEIGHT_MM = 1.0                # a p300 tip SEALS on a flat floor below ~0.8 mm
MIX_DEPTH_FRACTION = 0.5               # of the liquid left at the stroke's bottom
ASPIRATE_MIN_HEIGHT_MM = 0.8           # same floor problem, same fix
ASPIRATE_DEPTH_FRACTION = 0.5          # of the liquid left AFTER the stroke
ASPIRATE_MARGIN_MM = 0.3               # liquid that MUST remain over the tip at the end
                                       # of any aspirate. draw_at raises below this
                                       # rather than silently drawing air.
BLOWOUT_ABOVE_MM = 1.0                 # blow out above the surface, never under it
# --- droplet control -----------------------------------------------------------------
BLOW_OUT = True                        # clears the INSIDE of the tip
TOUCH_TIP = True                       # sheds the drop on the OUTSIDE
TOUCH_RADIUS = 0.8                     # fraction of well radius
TOUCH_V_OFFSET_MM = -1.5               # below the well top
# --- multi-dispense ------------------------------------------------------------------
# One aspirate serves several wells, with a touch-tip after every dispense. The default
# disposal volume is the pipette's own minimum (the Opentrons distribute() default): it
# rides on top of the load and is blown into the trash, so the last well out of the tip
# is as accurate as the first. build_layout() adds that waste to the reservoir totals.
# The 10 uL spikes are the exception - 2 x 10 already fills the p20 to its 20 uL brim,
# so they run with no disposal volume at all.


# =======================================================================================
# SECTION 4 - LAYOUT   (fixed blocks; nothing here needs editing)
# =======================================================================================
ARM_PLAIN, ARM_GUA = 'plain', 'guaiacol'         # the two NNBT arms
MIX_NNBT, MIX_GUA = 'nnbt', 'gua'                # which reservoir reagent a group takes
MIX_NO_NNBT, MIX_NO_NNBT_GUA = 'no_nnbt', 'no_nnbt_gua'

MIX_RESERVOIR = {MIX_NNBT: RES_NNBT, MIX_GUA: RES_GUAIACOL,
                 MIX_NO_NNBT: RES_NO_NNBT, MIX_NO_NNBT_GUA: RES_NO_NNBT_GUA}
MIX_LABEL = {MIX_NNBT: 'NNBT mix', MIX_GUA: 'guaiacol mix',
             MIX_NO_NNBT: 'no-NNBT mix', MIX_NO_NNBT_GUA: 'no-NNBT guaiacol mix'}


def _tube(index):
    """Tube-rack well for the index-th enzyme (0-based): 4 rows per column."""
    return f'{ROW_LETTERS[index % TUBERACK_ROWS]}{index // TUBERACK_ROWS + 1}'


def _dil(index):
    """Dilution-plate well for the index-th dilution (0-based), column-then-row."""
    return f'{ROW_LETTERS[index % WELLS_PER_COLUMN]}{index // WELLS_PER_COLUMN + 1}'


def _block_wells(block, slot):
    """The 3 side-by-side wells of group `slot` (0-7) in 3-column `block` (0-based)."""
    col0 = block * BLOCK_COLS + 1
    return [f'{ROW_LETTERS[slot]}{col0 + k}' for k in range(REPLICATES)]


def to_um(mg_ml, mw_kda):
    """Bradford mass concentration -> molarity.  uM = mg/mL / kDa x 1000."""
    if mw_kda <= 0:
        raise ValueError('molar mass must be > 0 kDa')
    return mg_ml / mw_kda * 1000.0


def dilution_recipe(stock_conc, target_conc, total=DILUTION_WELL_VOL_UL):
    """(stock_uL, buffer_uL) making `total` at `target_conc`. Raises rather than
    silently producing a dilution the pipettes cannot deliver."""
    if stock_conc <= 0:
        raise ValueError('stock concentration must be > 0')
    stock = round(total * target_conc / stock_conc, 1)
    if stock > total + 1e-9:
        raise ValueError(
            f'stock at {stock_conc:.3g} is WEAKER than the target {target_conc:.3g}: '
            f'it would need {stock:.1f} uL in a {total:g} uL well. Lower the target.')
    if stock < P20_MIN:
        raise ValueError(
            f'stock at {stock_conc:.3g} is so strong only {stock:.2f} uL is needed - '
            f'below the p20 minimum of {P20_MIN:g} uL. Pre-dilute it by hand.')
    return stock, round(total - stock, 1)


def build_layout(enzymes, heat_um, target_um=None, nc3_index=0):
    """The single source of truth: dilutions, both plate maps, reagent volumes.
    `enzymes` is [{'name','mg_ml','mw_kda'}]. Pure Python - runs on your laptop."""
    n = len(enzymes)
    if n > MAX_ENZYMES:
        raise ValueError(f'{n} enzymes but only {MAX_ENZYMES} fit: each one needs a '
                         'well in both the plain and the guaiacol block.')
    enz_um = [{'name': e['name'], 'conc': to_um(e['mg_ml'], e['mw_kda']),
               'mg_ml': e['mg_ml'], 'mw_kda': e['mw_kda']} for e in enzymes]

    # ---- target: the common molar concentration -------------------------------------
    # The heat-inactivated control is diluted like a sample, so it counts too.
    target = float(target_um) if target_um else (
        min([e['conc'] for e in enz_um] + [heat_um]) * ENZ_TARGET_AUTO_FRACTION)

    # ---- dilution plate: enzymes, NC3's enzyme, heat-inactivated, lac curve, ABTS ----
    dilutions, idx = [], 0

    def add_dil(name, kind, address, stock_conc, note=''):
        nonlocal idx
        stock, buf = dilution_recipe(stock_conc, target)
        d = {'name': name, 'well': _dil(idx), 'src_kind': kind, 'src': address,
             'stock_conc': stock_conc, 'stock_vol': stock, 'buffer_vol': buf,
             'note': note}
        dilutions.append(d)
        idx += 1
        return d

    enz_dils = [add_dil(e['name'], 'tube', _tube(i), e['conc'])
                for i, e in enumerate(enz_um)]
    nc3_dil = add_dil(f'NC3 ({enz_um[nc3_index]["name"]})', 'tube', _tube(nc3_index),
                      enz_um[nc3_index]['conc'], 'second aliquot, keeps NC3 off the '
                      'enzyme well so neither is over-drawn')
    heat_dil = add_dil('NC2 heat-inactivated', 'tube', TUBE_HEAT_INACT, heat_um)

    # lactaldehyde curve: strongest first, each point serial-diluted from the previous,
    # so a source is always finished before anything draws from it.
    spike_factor = REACTION_VOL_UL / SPIKE_VOL_UL          # 15x
    lac_dils, src_conc, src = [], float(LAC_STOCK_UM), ('tube', TUBE_LACTALDEHYDE)
    for um in sorted(LAC_GRADIENT_UM, reverse=True):
        want = um * spike_factor
        stock, buf = dilution_recipe(src_conc, want)
        d = {'name': f'Lac-std {um:g} uM', 'well': _dil(idx), 'src_kind': src[0],
             'src': src[1], 'stock_conc': src_conc, 'stock_vol': stock,
             'buffer_vol': buf, 'note': 'serial' if src[0] == 'dil' else ''}
        dilutions.append(d)
        lac_dils.append(d)
        if LAC_SERIAL:                                 # chain: next point comes from
            src_conc, src = want, ('dil', d['well'])   # this one. See LAC_SERIAL.
        idx += 1

    # ABTS standards: pre-made, just moved off the rack before it leaves the deck.
    abts_dils = []
    for conc, tube in zip(ABTS_STD_MG_ML, ABTS_STD_TUBES):
        d = {'name': f'ABTS-std {conc:g} mg/mL', 'well': _dil(idx), 'src_kind': 'tube',
             'src': tube, 'stock_conc': conc, 'stock_vol': 0.0, 'buffer_vol': 0.0,
             'note': 'pre-made, transferred not diluted'}
        dilutions.append(d)
        abts_dils.append(d)
        idx += 1

    # ---- NNBT plate: four blocks of eight groups ------------------------------------
    def ctrl_block(block, arm, mix, nomix):
        """One control block: NC1, NC2, NC3, then the 5 lactaldehyde standards."""
        groups = [
            {'label': 'NC1 buffer only', 'arm': arm, 'mix': mix,
             'src': ('res', RES_BUFFER)},
            {'label': heat_dil['name'], 'arm': arm, 'mix': mix,
             'src': ('dil', heat_dil['well'])},
            {'label': nc3_dil['name'], 'arm': arm, 'mix': nomix,
             'src': ('dil', nc3_dil['well'])},
        ] + [{'label': d['name'], 'arm': arm, 'mix': mix, 'src': ('dil', d['well'])}
             for d in lac_dils]
        for slot, g in enumerate(groups):
            g['wells'] = _block_wells(block, slot)
        return groups

    def enz_block(block, arm, mix):
        groups = [{'label': d['name'], 'arm': arm, 'mix': mix,
                   'src': ('dil', d['well'])} for d in enz_dils]
        for slot, g in enumerate(groups):
            g['wells'] = _block_wells(block, slot)
        return groups

    nnbt = (ctrl_block(0, ARM_PLAIN, MIX_NNBT, MIX_NO_NNBT)
            + ctrl_block(1, ARM_GUA, MIX_GUA, MIX_NO_NNBT_GUA)
            + enz_block(2, ARM_PLAIN, MIX_NNBT)
            + enz_block(3, ARM_GUA, MIX_GUA))

    # ---- ABTS plate: NC1, NC2, the 6 standards, then the enzymes --------------------
    abts = [{'label': 'NC1 buffer only', 'src': ('res', RES_BUFFER)},
            {'label': heat_dil['name'], 'src': ('dil', heat_dil['well'])}]
    abts += [{'label': d['name'], 'src': ('dil', d['well'])} for d in abts_dils]
    for slot, g in enumerate(abts):
        g['wells'] = _block_wells(0, slot)
    for slot, d in enumerate(enz_dils):
        abts.append({'label': d['name'], 'src': ('dil', d['well']),
                     'wells': _block_wells(1, slot)})
    abts_columns = max(int(w[1:]) for g in abts for w in g['wells'])

    # ---- reagent volumes ------------------------------------------------------------
    # Blocks 1-2 are filled one well at a time (NC3 needs its own mix), so they cost
    # only the wells actually used. Blocks 3-4 are 8-channel, so a part-filled column
    # still costs a full column.
    per_mix = {}
    ctrl_wells = {}
    for g in nnbt[:16]:                                   # the two control blocks
        per_mix[g['mix']] = per_mix.get(g['mix'], 0.0) + REACTION_MIX_VOL_UL * REPLICATES
        ctrl_wells[g['mix']] = ctrl_wells.get(g['mix'], 0) + REPLICATES
    enz_cols = BLOCK_COLS if n else 0
    per_mix[MIX_NNBT] = per_mix.get(MIX_NNBT, 0.0) + \
        REACTION_MIX_VOL_UL * WELLS_PER_COLUMN * enz_cols
    per_mix[MIX_GUA] = per_mix.get(MIX_GUA, 0.0) + \
        REACTION_MIX_VOL_UL * WELLS_PER_COLUMN * enz_cols

    # MULTI-DISPENSE WASTE. Every aspirate carries a disposal volume that is blown into
    # the trash, so the reservoir has to hold more than the wells actually consume.
    # Same batching arithmetic as multi_dispense() in run() - keep the two in step.
    for key, wells_n in ctrl_wells.items():               # one well at a time
        per_mix[key] += P300_MIN * p300_loads(wells_n, REACTION_MIX_VOL_UL)
    if enz_cols:                                          # 8-channel: waste x 8 rows
        for key in (MIX_NNBT, MIX_GUA):
            per_mix[key] += (P300_MIN * p300_loads(enz_cols, REACTION_MIX_VOL_UL)
                             * WELLS_PER_COLUMN)

    nnbt_columns = 2 * BLOCK_COLS + (2 * BLOCK_COLS if n else 0)
    purpald_total = (PURPALD_VOL_UL * WELLS_PER_COLUMN * nnbt_columns
                     + P300_MIN * p300_loads(nnbt_columns, PURPALD_VOL_UL)
                     * WELLS_PER_COLUMN)

    buffer_total = (sum(d['buffer_vol'] for d in dilutions)
                    + SPIKE_VOL_UL * REPLICATES * 2      # NC1 on both NNBT arms
                    + SPIKE_VOL_UL * REPLICATES)         # NC1 on the ABTS plate
    return {
        'target': target, 'enzymes': enz_um, 'dilutions': dilutions,
        'enz_dils': enz_dils, 'nc3_dil': nc3_dil, 'heat_dil': heat_dil,
        'lac_dils': lac_dils, 'abts_dils': abts_dils,
        'nnbt': nnbt, 'abts': abts,
        'nnbt_columns': nnbt_columns,
        'abts_columns': abts_columns,
        'per_mix': per_mix, 'buffer_total': buffer_total,
        'purpald_total': purpald_total,
        'abts_total': ABTS_MIX_VOL_UL * WELLS_PER_COLUMN * abts_columns,
    }


def p300_loads(n_dests, vol, disposal=P300_MIN):
    """How many p300 aspirates multi_dispense() needs for n_dests wells of `vol`."""
    if not n_dests:
        return 0
    per_load = max(1, int((P300_MAX - disposal) // vol))
    return math.ceil(n_dests / per_load)


def spread(total_ul, wells, reagent):
    """Split a reagent evenly over the fewest reservoir wells that hold it."""
    need = max(1, math.ceil(total_ul / RES_USABLE_UL))
    if need > len(wells):
        raise ValueError(f'{reagent} needs {total_ul / 1000:.1f} mL = {need} wells, '
                         f'but only {len(wells)} are configured.')
    return {w: total_ul / need for w in wells[:need]}


# --- preflight: refuse to load on a rule we have already been burned by ----------------
LABWARE_HEIGHT_MM = {'opentrons_96_tiprack_20ul': 64.7,
                     'opentrons_96_tiprack_300ul': 64.5, TUBERACK: 54.0,
                     RESERVOIR_LOADNAME: 31.4, NNBT_PLATE: 14.3, 'TRASH': 999.0}


def check_deck():
    """Slot 9 can never be a single-nozzle target (trash in 12 is north of it) and
    nothing tall may sit north of one. opentrons_simulate does NOT catch the trash
    case - it passed the layout the robot then refused."""
    occupied = {str(HS_SLOT): NNBT_PLATE, SLOT_SWAP: TUBERACK,
                SLOT_TIPRACK_20: 'opentrons_96_tiprack_20ul',
                SLOT_DILUTION: DILUTION_PLATE,
                SLOT_TIPRACK_300: 'opentrons_96_tiprack_300ul',
                SLOT_RESERVOIR: RESERVOIR_LOADNAME, '12': 'TRASH'}
    problems = []
    for slot, what in occupied.items():
        if slot == '12':
            continue
        north = str(int(slot) + 3)
        if north == '12':
            problems.append(f'{what} in slot {slot} has the fixed trash bin north of '
                            'it - this is the 2026-09-04 failure. Move it.')
        elif north in occupied and LABWARE_HEIGHT_MM.get(occupied[north], 0) >= 50:
            problems.append(f'{what} in slot {slot} has {occupied[north]} north of it '
                            '- too tall for a single-nozzle move.')
    if problems:
        raise ValueError('DECK RULE VIOLATION:\n  - ' + '\n  - '.join(problems))


def check_mix_geometry():
    """A tip above the liquid draws air; a tip on the floor seals. Both have happened."""
    area = math.pi * (6.85 / 2) ** 2
    stroke = min(DILUTION_WELL_VOL_UL * MIX_STROKE_FRACTION, P300_MAX)
    surface = (DILUTION_WELL_VOL_UL - stroke) / area
    z = max(MIX_MIN_HEIGHT_MM, surface * MIX_DEPTH_FRACTION)
    if z >= surface:
        raise ValueError(f'mix tip at {z:.2f} mm but a {stroke:g} uL stroke drops the '
                         f'surface to {surface:.2f} mm - it would draw AIR.')
    if z < 0.8:
        raise ValueError(f'mix tip at {z:.2f} mm above a FLAT floor - a p300 tip seals '
                         'below ~0.8 mm and moves nothing.')


check_deck()
check_mix_geometry()


# Colours for the Opentrons app's labware map. The app shows a coloured, named liquid
# in every tube, reservoir well and plate well, so the setup screen doubles as the
# "what goes where" sheet.
COLOUR = {
    'enzyme':    '#2f7d5f',            # green  - the enzymes under test
    'nc1':       '#78706f',            # grey   - buffer only
    'nc2':       '#3f6d8f',            # blue   - heat-inactivated
    'nc3':       '#7b4fa0',            # purple - active enzyme, no NNBT
    'lac':       '#c2621a',            # amber  - lactaldehyde standards
    'abts_std':  '#1b7f8f',            # teal   - ABTS protein standards
    'nnbt':      '#8e24aa',            # violet - NNBT reaction mix
    'gua':       '#b8860b',            # gold   - guaiacol reaction mix
    'no_nnbt':   '#c2185b',            # pink   - reaction mix without NNBT
    'buffer':    '#0288d1',            # cyan   - assay buffer
    'purpald':   '#c62828',            # red    - Purpald
    'abts_mix':  '#00695c',            # dark teal
    'water':     '#90a4ae',            # pale   - milliQ
    'matrix':    '#5d4037',            # brown  - MALDI matrix
}


def colour_of(name):
    """Pick a colour from a group/dilution label."""
    if name.startswith('NC1'):
        return COLOUR['nc1']
    if name.startswith('NC2'):
        return COLOUR['nc2']
    if name.startswith('NC3'):
        return COLOUR['nc3']
    if name.startswith('Lac-std'):
        return COLOUR['lac']
    if name.startswith('ABTS-std'):
        return COLOUR['abts_std']
    return COLOUR['enzyme']


# =======================================================================================
# SECTION 5 - REPORTING   (run log and `python3 <this file>`)
# =======================================================================================
SYMBOL = {MIX_NNBT: 'N', MIX_GUA: 'G', MIX_NO_NNBT: 'x', MIX_NO_NNBT_GUA: 'y'}


def render(layout, n_enz, maldi_on=False, interval=30, maldi_row='A'):
    """Everything you need to set the deck up and read the plates afterwards."""
    out = ['=' * 78,
           f'  {n_enz} enzymes x 2 arms (plain + guaiacol), 3 negative controls, '
           f'{len(LAC_GRADIENT_UM)} lactaldehyde standards',
           f'  every enzyme diluted to {layout["target"]:.2f} uM',
           '=' * 78, '']

    # --- NNBT plate map ---
    pos = {w: i + 1 for i, g in enumerate(layout['nnbt']) for w in g['wells']}
    mix = {w: g['mix'] for g in layout['nnbt'] for w in g['wells']}
    out.append('  NNBT PLATE (slot 1)   N=NNBT  G=guaiacol  x=no-NNBT  y=no-NNBT+gua')
    out.append('       ' + ''.join(f'{c:>5}' for c in range(1, 13)))
    for r in ROW_LETTERS:
        cells = [f'{SYMBOL[mix[f"{r}{c}"]]}{pos[f"{r}{c}"]:02d}'.rjust(5)
                 if f'{r}{c}' in pos else '    .' for c in range(1, 13)]
        out.append(f'  {r}  ' + ''.join(cells))
    out.append('')
    out.append('   #  wells          content')
    for i, g in enumerate(layout['nnbt']):
        out.append(f'  {i + 1:02d}  {",".join(g["wells"]):<14} {g["label"]} '
                   f'[{MIX_LABEL[g["mix"]]}]')

    # --- ABTS plate map ---
    out += ['', '-' * 78, f'  ABTS PLATE (slot {SLOT_SWAP}, swapped in after the dilutions)',
            '   #  wells          content']
    for i, g in enumerate(layout['abts']):
        out.append(f'  {i + 1:02d}  {",".join(g["wells"]):<14} {g["label"]}')

    # --- dilution plate ---
    out += ['', '-' * 78,
            f'  DILUTION PLATE (slot 6), {DILUTION_WELL_VOL_UL:g} uL per well',
            '  well  from        stock uL  buffer uL  content']
    for d in layout['dilutions']:
        origin = ('tube ' if d['src_kind'] == 'tube' else 'dil  ') + d['src']
        out.append(f'  {d["well"]:<5} {origin:<11} {d["stock_vol"]:>8.1f} '
                   f'{d["buffer_vol"]:>10.1f}  {d["name"]}'
                   f'{"  [" + d["note"] + "]" if d["note"] else ""}')

    # --- what to pour ---
    needs = {RES_BUFFER: ('assay buffer', layout['buffer_total'])}
    for key, total in layout['per_mix'].items():
        for w, v in spread(total, MIX_RESERVOIR[key], MIX_LABEL[key]).items():
            needs[w] = (MIX_LABEL[key], v)
    for w, v in spread(layout['purpald_total'], RES_PURPALD, 'Purpald').items():
        needs[w] = ('Purpald reagent', v)
    for w, v in spread(layout['abts_total'], RES_ABTS, 'ABTS mix').items():
        needs[w] = ('ABTS reaction mix', v)
    if maldi_on:
        rounds = int(INCUBATION_MIN // interval) + 1
        spots = 2 * (3 + n_enz) * rounds
        needs[RES_WATER] = ('milliQ (MALDI 1:5)', MALDI_WATER_UL * spots)
        needs[RES_MATRIX] = ('MALDI matrix', MALDI_MATRIX_UL * spots)
    out += ['', '-' * 78, f'  RESERVOIR (slot {SLOT_RESERVOIR}) - pour these, '
            f'{RES_DEAD_UL / 1000:g} mL dead volume already included']
    for w in sorted(needs, key=lambda x: int(x[1:])):
        name, vol = needs[w]
        out.append(f'    {w:<4} {name:<24} {(vol + RES_DEAD_UL) / 1000:>5.1f} mL')

    # --- tubes ---
    out += ['', f'  TUBE RACK (slot {SLOT_SWAP}, removed after the dilutions)']
    tubes = {}                                    # NC3 shares a tube with its enzyme
    for d in layout['dilutions']:
        if d['src_kind'] != 'tube':
            continue
        vol, names = tubes.get(d['src'], (0.0, []))
        tubes[d['src']] = (vol + (d['stock_vol'] or 40.0), names + [d['name']])
    for tube in sorted(tubes, key=lambda t: (int(t[1:]), t[0])):
        vol, names = tubes[tube]
        out.append(f'    {tube:<4} {" + ".join(names):<34} >= {vol + 300:.0f} uL')

    # --- MALDI ---
    if maldi_on:
        rounds = int(INCUBATION_MIN // interval) + 1
        start = MALDI_ROWS.index(maldi_row)
        out += ['', '-' * 78,
                f'  MALDI - {rounds} rounds every {interval} min, one spot per condition',
                '  round   t (min)   plain row   guaiacol row']
        for r in range(rounds):
            if start + 2 * r + 1 >= len(MALDI_ROWS):
                raise ValueError(f'MALDI needs {2 * rounds} rows from {maldi_row}; the '
                                 f'target has {len(MALDI_ROWS)}. Start higher up or '
                                 'use a longer interval.')
            out.append(f'  {r + 1:<7} {r * interval:>7}   '
                       f'{MALDI_ROWS[start + 2 * r]:<11} {MALDI_ROWS[start + 2 * r + 1]}')
        out.append(f'  columns 1..{3 + n_enz}: NC1, NC2, NC3, then the enzymes')
    return out


# =======================================================================================
# SECTION 6 - THE PROTOCOL
# =======================================================================================
def run(protocol):
    prm = protocol.params
    n = prm.n_enz
    rounds = int(INCUBATION_MIN // prm.maldi_interval) + 1 if prm.maldi_on else 0
    # rounds = 1

    # ---- the batch, straight from the app ------------------------------------------
    batch = [{'name': ENZYME_NAMES[i], 'mg_ml': getattr(prm, f'mg_ml_{i + 1}'),
              'mw_kda': getattr(prm, f'mw_{i + 1}')} for i in range(n)]
    heat_um = to_um(prm.heat_mg_ml, prm.heat_mw)
    layout = build_layout(batch, heat_um, prm.target_um, prm.nc3_enz - 1)

    for e in layout['enzymes']:                       # the conversion, on the record
        protocol.comment(f'  {e["name"]}: {e["mg_ml"]:g} mg/mL / {e["mw_kda"]:g} kDa '
                         f'= {e["conc"]:.2f} uM')
    for line in render(layout, n, prm.maldi_on, prm.maldi_interval, prm.maldi_row):
        protocol.comment(line)

    # ---- deck ----------------------------------------------------------------------
    hs = protocol.load_module('heaterShakerModuleV1', HS_SLOT)
    nnbt_plate = hs.load_adapter(HS_ADAPTER).load_labware(NNBT_PLATE)
    tuberack = protocol.load_labware(TUBERACK, SLOT_SWAP)
    dil_plate = protocol.load_labware(DILUTION_PLATE, SLOT_DILUTION)
    reservoir = protocol.load_labware(RESERVOIR_LOADNAME, SLOT_RESERVOIR)
    tr20 = protocol.load_labware('opentrons_96_tiprack_20ul', SLOT_TIPRACK_20)
    tr300 = protocol.load_labware('opentrons_96_tiprack_300ul', SLOT_TIPRACK_300)
    abts_plate = protocol.load_labware(ABTS_PLATE, protocol_api.OFF_DECK)
    maldi_target = (protocol.load_labware(MALDI_PLATE, protocol_api.OFF_DECK)
                    if prm.maldi_on else None)
    maldi_dil = (protocol.load_labware(DILUTION_PLATE, protocol_api.OFF_DECK)
                 if prm.maldi_on else None)

    p20 = protocol.load_instrument('p20_multi_gen2', 'right', tip_racks=[tr20])
    p300 = protocol.load_instrument('p300_multi_gen2', 'left', tip_racks=[tr300])
    p20.configure_nozzle_layout(style=SINGLE, start='H1', tip_racks=[tr20])
    p300.configure_nozzle_layout(style=SINGLE, start='H1', tip_racks=[tr300])
    p20.flow_rate.aspirate = p20.flow_rate.dispense = P20_FLOW_UL_S    # how fast liquid
    p300.flow_rate.aspirate = p300.flow_rate.dispense = P300_FLOW_UL_S # moves in/out
    p20.flow_rate.blow_out = P20_BLOWOUT_UL_S       # kept gentle on purpose
    p300.flow_rate.blow_out = P300_BLOWOUT_UL_S
    hs.close_labware_latch()

    # ---- LIQUID MAP -----------------------------------------------------------------
    # define_liquid + load_liquid makes the app draw a coloured, named liquid in every
    # tube, reservoir well and plate well. The setup screen then tells you exactly what
    # to put where, and how much, without reading this file.
    liq = {}                                       # name -> Liquid, defined once

    def liquid(name, colour, note=''):
        if name not in liq:
            liq[name] = protocol.define_liquid(name, note or name, colour)
        return liq[name]

    # tubes: this is the "what goes in each tube" answer
    tube_needs = {}
    for d in layout['dilutions']:
        if d['src_kind'] == 'tube':
            vol, names = tube_needs.get(d['src'], (0.0, []))
            tube_needs[d['src']] = (vol + (d['stock_vol'] or 40.0), names + [d['name']])
    for tube, (vol, names) in tube_needs.items():
        label = ' + '.join(dict.fromkeys(names))
        tuberack[tube].load_liquid(liquid(label, colour_of(names[0])), vol + 300)

    # reservoir: every reagent, with the volume to pour
    res_liquids = [(RES_BUFFER, 'Assay buffer', 'buffer', layout['buffer_total'])]
    for key, total in layout['per_mix'].items():
        tag = 'nnbt' if key == MIX_NNBT else 'gua' if key == MIX_GUA else 'no_nnbt'
        for w, v in spread(total, MIX_RESERVOIR[key], MIX_LABEL[key]).items():
            res_liquids.append((w, MIX_LABEL[key].capitalize(), tag, v))
    for w, v in spread(layout['purpald_total'], RES_PURPALD, 'Purpald').items():
        res_liquids.append((w, 'Purpald reagent', 'purpald', v))
    for w, v in spread(layout['abts_total'], RES_ABTS, 'ABTS mix').items():
        res_liquids.append((w, 'ABTS reaction mix', 'abts_mix', v))
    if prm.maldi_on:
        spots = 2 * (3 + n) * rounds
        res_liquids.append((RES_WATER, 'milliQ (MALDI 1:5)', 'water',
                            MALDI_WATER_UL * spots))
        res_liquids.append((RES_MATRIX, 'MALDI matrix', 'matrix',
                            MALDI_MATRIX_UL * spots))
    for well, name, tag, vol in res_liquids:
        reservoir[well].load_liquid(liquid(name, COLOUR[tag]), vol + RES_DEAD_UL)

    # plates: start empty, so volume 0 - the colour and name are the point
    for d in layout['dilutions']:
        dil_plate[d['well']].load_liquid(liquid(d['name'], colour_of(d['name'])), 0)
    for g in layout['nnbt']:
        tag = ' + guaiacol' if g['arm'] == ARM_GUA else ''
        for w in g['wells']:
            nnbt_plate[w].load_liquid(
                liquid(g['label'] + tag, colour_of(g['label'])), 0)
    for g in layout['abts']:
        for w in g['wells']:
            abts_plate[w].load_liquid(
                liquid('ABTS: ' + g['label'], colour_of(g['label'])), 0)

    # ---- tips: one running index per rack, so single and column picks never collide --
    racks = {'20': tr20.wells(), '300': tr300.wells()}
    start = {'20': racks['20'].index(tr20.well(f'{prm.tip20_row}{prm.tip20_col}')),
             '300': racks['300'].index(tr300.well(f'A{prm.tip300_col}'))}
    nxt = dict(start)

    def take(key, pip, slot, count=1):
        if nxt[key] + count > len(racks[key]):        # ran out - ask for a fresh rack
            protocol.pause(f'{key} uL tips used up. Put a FULL rack in slot {slot} '
                           'and resume.')
            nxt[key] = 0
        well = racks[key][nxt[key]]                   # row A when count == 8
        nxt[key] += count
        pip.pick_up_tip(well)

    def slow():
        """Gentle p20, for placing 1 uL droplets on the MALDI target."""
        p20.flow_rate.aspirate = p20.flow_rate.dispense = MALDI_FLOW_UL_S

    def fast():
        p20.flow_rate.aspirate = p20.flow_rate.dispense = P20_FLOW_UL_S

    def pick20():
        take('20', p20, SLOT_TIPRACK_20)

    def pick300():
        take('300', p300, SLOT_TIPRACK_300)

    def pick300_column():
        while nxt['300'] % WELLS_PER_COLUMN:          # skip a part-used column
            nxt['300'] += 1
        take('300', p300, SLOT_TIPRACK_300, WELLS_PER_COLUMN)

    # ---- tip budget: warn BEFORE the run rather than stalling mid-incubation --------
    per_arm = 3 + n                                   # 3 NCs + the enzymes
    maldi_tips = rounds * 2 * (per_arm + math.ceil(per_arm / MALDI_MATRIX_BATCH))
    tips20 = (len(layout['dilutions']) + len(layout['nnbt']) + len(layout['abts'])
              + maldi_tips + 2)
    have20 = len(racks['20']) - start['20']
    protocol.comment(f'  20 uL tips needed ~{tips20}, available from '
                     f'{prm.tip20_row}{prm.tip20_col}: {have20}')
    if tips20 > have20:
        protocol.comment(f'  *** NOT ENOUGH 20 uL TIPS: have {have20}, need ~{tips20}. '
                         f'The run will PAUSE for a fresh rack {math.ceil((tips20 - have20) / 96)} '
                         'time(s) - have them open and ready. ***')

    # ---- MALDI fits? check now, not two hours into the incubation -------------------
    if prm.maldi_on:
        per_round = 2 * (3 + n)                        # both arms, 3 NCs + n enzymes
        if per_round * rounds > 96:
            raise ValueError(
                f'MALDI needs {per_round} dilution wells per round x {rounds} rounds = '
                f'{per_round * rounds}, but the plate has 96. Use a longer interval '
                f'(>= {math.ceil(INCUBATION_MIN / (96 / per_round - 1)):.0f} min) or '
                'fewer enzymes.')
        if MALDI_ROWS.index(prm.maldi_row) + 2 * rounds > len(MALDI_ROWS):
            raise ValueError(
                f'MALDI needs {2 * rounds} rows from {prm.maldi_row}; the target has '
                f'{len(MALDI_ROWS)}. Start higher up or use a longer interval.')

    # ---- small helpers --------------------------------------------------------------
    def area(well):
        d = getattr(well, 'diameter', None)
        return math.pi * (d / 2) ** 2 if d else well.length * well.width

    def height(well, vol):
        return vol / area(well)

    def mix_at(well, full, stroke):
        """Submerged for the WHOLE stroke: the level that matters is the one at the
        BOTTOM of the stroke, not at rest."""
        low = height(well, max(0.0, full - stroke))
        return well.bottom(z=min(max(MIX_MIN_HEIGHT_MM, low * MIX_DEPTH_FRACTION),
                                 well.depth - 1.0))

    def draw_at(well, left, stroke):
        """Submerged for the whole aspirate: the level AFTER it, not before."""
        after = max(0.0, left - stroke)
        surface = height(well, after)
        z = min(max(ASPIRATE_MIN_HEIGHT_MM, surface * ASPIRATE_DEPTH_FRACTION),
                well.depth - 1.0)
        # The half-depth rule keeps the tip under the surface, but ASPIRATE_MIN_HEIGHT_MM
        # overrides it once the well is nearly empty - and then the tip can end up AT or
        # ABOVE the liquid. Refuse rather than quietly aspirate air. This fires during a
        # laptop simulate, before the robot touches anything.
        if z >= surface - ASPIRATE_MARGIN_MM:
            raise ValueError(
                f'{well}: drawing {stroke:.1f} uL from {left:.1f} uL leaves the surface '
                f'at {surface:.2f} mm but the tip sits at {z:.2f} mm (needs '
                f'{ASPIRATE_MARGIN_MM:g} mm clearance) - it would aspirate AIR. '
                f'Raise DILUTION_WELL_VOL_UL, or split this transfer.')
        return well.bottom(z=z)

    def as_well(t):
        return t.labware.as_well() if hasattr(t, 'labware') else t

    def touch(pip, target):
        w = as_well(target)
        if TOUCH_TIP and w is not None:
            pip.touch_tip(w, radius=TOUCH_RADIUS, v_offset=TOUCH_V_OFFSET_MM)

    def clear(pip, target, over=None):
        """Blow out above any liquid, then shed the outside drop."""
        w = as_well(target)
        if w is None:
            return
        if BLOW_OUT:
            z = (height(w, over) + BLOWOUT_ABOVE_MM) if over else BLOWOUT_ABOVE_MM + 1.0
            pip.blow_out(w.bottom(z=min(z, w.depth - 1.0)))
        touch(pip, w)

    def transfer(vol, src, dest, pip=None, mix_well=None, new_tip=True, keep=False):
        """One single-nozzle transfer, split into strokes if bigger than the pipette."""
        pip = pip or (p20 if vol <= P20_MAX else p300)
        if new_tip:
            pick20() if pip is p20 else pick300()
        cap = P20_MAX if pip is p20 else P300_MAX
        left, done = vol, 0.0
        while left > 1e-6:
            stroke = min(left, cap)
            pip.aspirate(stroke, src)
            touch(pip, src)                            # drop falls back into the source
            pip.dispense(stroke, dest)
            done, left = done + stroke, left - stroke
            if left > 1e-6:
                clear(pip, dest, over=done)            # each stroke lands in full
        if mix_well is not None:
            stroke = min(DILUTION_WELL_VOL_UL * MIX_STROKE_FRACTION, cap)
            pip.mix(MIX_REPS, stroke, mix_at(mix_well, DILUTION_WELL_VOL_UL, stroke))
            clear(pip, mix_well, over=DILUTION_WELL_VOL_UL)
        else:
            clear(pip, dest, over=done)
        if new_tip and not keep:
            pip.drop_tip()

    def multi_dispense(pip, vol, src, dests, height=REAGENT_HEIGHT_MM,
                       touch_height=None, disposal=None, touch_src=False):
        """MULTI-DISPENSE: one aspirate serves several wells, touch-tip after each.

        src         a well/location, or a callable given the load volume and returning
                    one (a dilution well whose surface is falling as we draw from it)
        disposal    uL drawn on top of the load and blown into the TRASH at the end of
                    each load, so the last well out of the tip is as accurate as the
                    first. Defaults to the pipette's own minimum - the Opentrons
                    distribute() default. Never returned to the reservoir: the tip has
                    touched destination wells by then.
        disposal=0  exact load, nothing to spare. The LAST dispense of each load then
                    gets the blow-out (there is nothing left behind to push out); the
                    others are placed at `touch_height` and shed by the touch-tip alone.
        """
        cap = P20_MAX if pip is p20 else P300_MAX
        disposal = (P20_MIN if pip is p20 else P300_MIN) if disposal is None else disposal
        touch_height = height if touch_height is None else touch_height
        per_load = max(1, int((cap - disposal + 1e-9) // vol))
        i = 0
        while i < len(dests):
            batch = dests[i:i + per_load]
            load = vol * len(batch) + disposal
            where = src(load) if callable(src) else src
            pip.aspirate(load, where)
            if touch_src:
                touch(pip, where)                  # drop falls back into the source
            for k, d in enumerate(batch):
                last = (k == len(batch) - 1)
                if disposal == 0 and last:
                    pip.dispense(vol, d.bottom(z=height))
                    clear(pip, d)                  # tip is empty: blow out, then touch
                else:
                    pip.dispense(vol, d.bottom(z=touch_height))
                    touch(pip, d)
            if disposal:
                pip.blow_out(protocol.fixed_trash)
            i += len(batch)

    def dil_well(name):
        return dil_plate[name].bottom(z=DILUTION_HEIGHT_MM)

    def source_of(spec, left, stroke=SPIKE_VOL_UL):
        """Where a spike load comes from: a bottomless reservoir, or a dilution well
        whose surface we follow down. `stroke` is the whole aspirate, which under
        multi-dispense is more than one well's worth."""
        kind, addr = spec
        if kind == 'res':
            return reservoir[addr]
        return draw_at(dil_plate[addr], left, stroke)

    def spike(groups, plate, taken):
        """SPIKE_VOL_UL into each replicate well, MULTI-DISPENSED, one tip per group.
        The p20 tops out at 20 uL, so one aspirate serves TWO wells and the third is
        its own stroke - 2 aspirates per triplicate instead of 3. No disposal volume
        fits in 2 x 10, so the second well of each pair takes the blow-out and the
        first is placed low and shed by its touch-tip.
        `taken` = how much has already been drawn from each dilution well."""
        for g in groups:
            pick20()
            state = {'left': DILUTION_WELL_VOL_UL - taken.get(g['src'][1], 0.0)}

            def src(load, g=g, state=state):
                where = source_of(g['src'], state['left'], load)
                state['left'] -= load
                return where

            multi_dispense(p20, SPIKE_VOL_UL, src, [plate[w] for w in g['wells']],
                           height=SPIKE_HEIGHT_MM, touch_height=SPIKE_LOW_HEIGHT_MM,
                           disposal=0.0, touch_src=True)
            if g['src'][0] == 'dil':
                taken[g['src'][1]] = taken.get(g['src'][1], 0.0) + SPIKE_VOL_UL * REPLICATES
            p20.drop_tip()

    # ===================================================================================
    # 1  buffer into every dilution well - empty wells, so one tip per pipette serves
    # ===================================================================================
    tipped = set()
    for d in layout['dilutions']:
        if d['buffer_vol'] <= 0:
            continue
        pip = p20 if d['buffer_vol'] <= P20_MAX else p300
        if pip not in tipped:
            pick20() if pip is p20 else pick300()
            tipped.add(pip)
        transfer(d['buffer_vol'], reservoir[RES_BUFFER], dil_well(d['well']),
                 pip=pip, new_tip=False)
    for pip in tipped:
        pip.drop_tip()

    # ===================================================================================
    # 2/3/4  stocks in and mixed; lactaldehyde curve; ABTS standards off the rack
    # ===================================================================================
    for d in layout['dilutions']:
        if d['stock_vol'] <= 0:                        # ABTS standards, handled below
            continue
        src = (tuberack[d['src']] if d['src_kind'] == 'tube'
               else draw_at(dil_plate[d['src']], DILUTION_WELL_VOL_UL, d['stock_vol']))
        transfer(d['stock_vol'], src, dil_well(d['well']),
                 mix_well=dil_plate[d['well']])

    for d in layout['abts_dils']:
        transfer(SPIKE_VOL_UL * REPLICATES + 10.0, tuberack[d['src']],
                 dil_well(d['well']))

    # ===================================================================================
    # 5  spikes into the NNBT plate - all four blocks
    # ===================================================================================
    # `taken` starts at what the DILUTION step already removed, not at zero. A lac
    # standard that fed the next point serially is NOT still at DILUTION_WELL_VOL_UL:
    # the strongest one gave away 3/4 of itself. Seeding this is what keeps the spike
    # tip following the real surface down instead of dipping above it.
    taken = {}
    for d in layout['dilutions']:
        if d['src_kind'] == 'dil' and d['stock_vol'] > 0:
            taken[d['src']] = taken.get(d['src'], 0.0) + d['stock_vol']
    spike(layout['nnbt'], nnbt_plate, taken)

    # ===================================================================================
    # 6/7  tube rack out, ABTS plate in, spike it
    # ===================================================================================
    protocol.move_labware(tuberack, protocol_api.OFF_DECK, use_gripper=False)
    protocol.move_labware(abts_plate, SLOT_SWAP, use_gripper=False)
    spike(layout['abts'], abts_plate, taken)

    # ===================================================================================
    # 8  reaction mixes into the NNBT plate
    #      control blocks one well at a time - NC3 takes a different mix from its
    #      neighbours and an 8-channel add cannot tell rows apart
    #      enzyme blocks 8-channel by column - all one mix
    # ===================================================================================
    res_of = {}                                        # which reservoir well per mix
    for key, total in layout['per_mix'].items():
        res_of[key] = list(spread(total, MIX_RESERVOIR[key], MIX_LABEL[key]))[0]

    slow_a = p300.flow_rate.aspirate * REACTION_MIX_FLOW_SCALE   # NNBT drips at speed
    slow_d = p300.flow_rate.dispense * REACTION_MIX_FLOW_SCALE
    fast_a, fast_d = p300.flow_rate.aspirate, p300.flow_rate.dispense
    p300.flow_rate.aspirate, p300.flow_rate.dispense = slow_a, slow_d

    # ONE FRESH TIP PER MIX. A tip that has dispensed NNBT mix carries NNBT residue,
    # which would put NNBT into NC3 and destroy the only control that tests for it.
    # The no-NNBT mixes go first, on a brand-new tip, before any NNBT is in play.
    order = [MIX_NO_NNBT, MIX_NO_NNBT_GUA, MIX_NNBT, MIX_GUA]
    for key in order:
        wells = [w for g in layout['nnbt'][:16] if g['mix'] == key for w in g['wells']]
        if not wells:
            continue
        pick300()                                      # fresh tip for this mix
        multi_dispense(p300, REACTION_MIX_VOL_UL, reservoir[res_of[key]],
                       [nnbt_plate[w] for w in wells])   # 2 wells per 300 uL aspirate
        p300.drop_tip()

    p300.configure_nozzle_layout(style=ALL, tip_racks=[tr300])   # back to 8-channel
    p300.flow_rate.aspirate, p300.flow_rate.dispense = slow_a, slow_d
    if n:                                              # enzyme blocks, whole columns
        for block, key in ((2, MIX_NNBT), (3, MIX_GUA)):
            pick300_column()                           # fresh column of tips per mix
            multi_dispense(p300, REACTION_MIX_VOL_UL, reservoir[res_of[key]],
                           [nnbt_plate[f'A{c}'] for c in
                            range(block * BLOCK_COLS + 1,
                                  block * BLOCK_COLS + 1 + BLOCK_COLS)])
            p300.drop_tip()
    p300.flow_rate.aspirate, p300.flow_rate.dispense = fast_a, fast_d

    # ===================================================================================
    # 9  ABTS mix - LAST liquid step: ABTS is kinetic and starts on contact
    # ===================================================================================
    abts_res = list(spread(layout['abts_total'], RES_ABTS, 'ABTS mix'))
    pick300_column()
    for i, c in enumerate(range(1, layout['abts_columns'] + 1)):
        src = abts_res[min(i // BLOCK_COLS, len(abts_res) - 1)]
        p300.aspirate(ABTS_MIX_VOL_UL, reservoir[src])
        p300.dispense(ABTS_MIX_VOL_UL, abts_plate[f'A{c}'].bottom(z=REAGENT_HEIGHT_MM))
        p300.touch_tip()
    p300.drop_tip()

    # ===================================================================================
    # 10  hand over: ABTS plate to the reader, seal the NNBT plate, MALDI labware in
    # ===================================================================================
    protocol.pause(f'TAKE THE ABTS PLATE from slot {SLOT_SWAP} to the reader NOW '
                   '(A414/A734) - it is already reacting. Then resume.')
    protocol.move_labware(abts_plate, protocol_api.OFF_DECK, use_gripper=False)
    if prm.maldi_on:
        protocol.move_labware(dil_plate, protocol_api.OFF_DECK, use_gripper=False)
        protocol.move_labware(maldi_dil, SLOT_DILUTION, use_gripper=False)
        protocol.move_labware(maldi_target, SLOT_SWAP, use_gripper=False)
    protocol.pause(f'Seal the NNBT plate for the {INCUBATION_MIN} min incubation at '
                   f'{INCUBATION_TEMP_C} degC, then resume.')

    # ===================================================================================
    # 11  incubate, pausing every interval to spot the MALDI target
    # ===================================================================================
    hs.set_and_wait_for_temperature(INCUBATION_TEMP_C)
    if not prm.maldi_on:
        hs.set_and_wait_for_shake_speed(INCUBATION_RPM)
        protocol.delay(minutes=INCUBATION_MIN)
        hs.deactivate_shaker()
    else:
        # One MALDI condition = one NNBT well, diluted 1:5 in milliQ, one spot. The
        # plain arm goes on one row, the guaiacol arm on the next.
        start_row = MALDI_ROWS.index(prm.maldi_row)
        conds = [(g['label'], g['arm'], g['wells'][0])          # sample replicate 1 only
                 for g in layout['nnbt']
                 if g['label'].startswith('NC') or g['src'][1] in
                 [d['well'] for d in layout['enz_dils']]]
        plain = [c for c in conds if c[1] == ARM_PLAIN]
        gua = [c for c in conds if c[1] == ARM_GUA]
        dil_i = [0]                               # running index into the fresh plate

        def next_maldi_well():
            if dil_i[0] >= 96:
                raise ValueError('the MALDI dilution plate is full: '
                                 f'{2 * len(plain)} wells per round x {rounds} rounds '
                                 '> 96. Use a longer interval or fewer enzymes.')
            w = _dil(dil_i[0])
            dil_i[0] += 1
            return w

        for r in range(rounds):
            if r:
                hs.set_and_wait_for_shake_speed(INCUBATION_RPM)
                protocol.delay(minutes=prm.maldi_interval)
                hs.deactivate_shaker()                 # never pipette a moving plate
            protocol.pause(f'MALDI round {r + 1}/{rounds} (t = {r * prm.maldi_interval} '
                           'min): UNSEAL the NNBT plate, then resume.')
            for arm_i, arm in enumerate((plain, gua)):
                row = MALDI_ROWS[start_row + 2 * r + arm_i]
                # NO blow-out, NO touch-tip, NO air gap on the target anywhere below:
                # air pushed through a 1 uL droplet sprays it onto neighbouring spots
                # and leaves a bubble, which is a hole in the crystal layer.
                for first in range(0, len(arm), MALDI_MATRIX_BATCH):
                    batch = arm[first:first + MALDI_MATRIX_BATCH]
                    # matrix pass: ONE tip, and it never touches a sample - otherwise
                    # it would carry sample back into the shared matrix well.
                    pick20()
                    slow()                            # droplet work: gentle or it splashes
                    for j in range(len(batch)):
                        spot = maldi_target[f'{row}{first + j + 1}']
                        p20.aspirate(MALDI_MATRIX_UL, reservoir[RES_MATRIX])
                        p20.dispense(MALDI_MATRIX_UL,
                                     spot.bottom(z=MALDI_SPOT_HEIGHT_MM))
                    fast()
                    p20.drop_tip()
                    # sample pass: fresh tip per condition, straight onto its matrix
                    for j, (label, _, src_well) in enumerate(batch):
                        col = first + j + 1
                        well = maldi_dil[next_maldi_well()]   # fresh dilution well
                        spot = maldi_target[f'{row}{col}']
                        pick20()
                        p20.aspirate(MALDI_WATER_UL, reservoir[RES_WATER])
                        p20.dispense(MALDI_WATER_UL, well.bottom(z=DILUTION_HEIGHT_MM))
                        p20.aspirate(MALDI_SAMPLE_UL,             # 5 uL of sample
                                     nnbt_plate[src_well].bottom(z=MALDI_DRAW_HEIGHT_MM))
                        p20.dispense(MALDI_SAMPLE_UL, well.bottom(z=DILUTION_HEIGHT_MM))
                        p20.mix(5, 10.0,
                                mix_at(well, MALDI_WATER_UL + MALDI_SAMPLE_UL, 10.0))
                        slow()                        # from here the tip is on steel
                        p20.aspirate(MALDI_SPOT_UL + 1.0,         # surplus stays in tip
                                     mix_at(well, MALDI_WATER_UL + MALDI_SAMPLE_UL, 2.0))
                        p20.dispense(MALDI_SPOT_UL,
                                     spot.bottom(z=MALDI_SPOT_HEIGHT_MM))
                        p20.mix(MALDI_MIX_REPS, MALDI_MIX_UL,     # mix in place, 2 uL
                                spot.bottom(z=MALDI_SPOT_HEIGHT_MM))
                        fast()
                        p20.drop_tip()
            protocol.pause('RESEAL the NNBT plate, then resume.')
        leftover = INCUBATION_MIN - (rounds - 1) * prm.maldi_interval
        if leftover > 0:
            hs.set_and_wait_for_shake_speed(INCUBATION_RPM)
            protocol.delay(minutes=leftover)
            hs.deactivate_shaker()
        protocol.pause(f'Take the MALDI target from slot {SLOT_SWAP}, dry it and apply '
                       'matrix. Spot key is in this log. Resume to finish.')

    # ===================================================================================
    # 12  Purpald, develop, hand off
    # ===================================================================================
    protocol.pause('Remove the seal from the NNBT plate, then resume for Purpald.')
    pur = list(spread(layout['purpald_total'], RES_PURPALD, 'Purpald'))[0]
    pick300_column()
    multi_dispense(p300, PURPALD_VOL_UL, reservoir[pur],          # 5 columns per load
                   [nnbt_plate[f'A{c}'] for c in range(1, layout['nnbt_columns'] + 1)])
    p300.drop_tip()

    hs.set_and_wait_for_temperature(DEVELOP_TEMP_C)
    hs.set_and_wait_for_shake_speed(DEVELOP_RPM)
    protocol.delay(minutes=DEVELOP_MIN)
    hs.deactivate_shaker()
    hs.deactivate_heater()
    hs.open_labware_latch()
    protocol.comment('Done - read A530. Plate maps are at the top of this log.')


# =======================================================================================
# SECTION 7 - LOCAL PREVIEW   `python3 <this file>` - no robot, no opentrons needed
# =======================================================================================
if __name__ == '__main__':
    demo = [{'name': 'Lac-01', 'mg_ml': 2.40, 'mw_kda': 65.0},
            {'name': 'Lac-02', 'mg_ml': 1.85, 'mw_kda': 65.0},
            {'name': 'Lac-03', 'mg_ml': 3.10, 'mw_kda': 70.0}]
    print('\n'.join(render(build_layout(demo, to_um(2.40, 65.0)),
                           len(demo), True, 30, 'A')))
