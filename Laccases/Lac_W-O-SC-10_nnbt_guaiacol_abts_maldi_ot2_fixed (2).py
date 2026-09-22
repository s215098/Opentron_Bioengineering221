"""
W-O-SC-10   NNBT/Purpald + guaiacol + ABTS + MALDI       Kristine Toft Johansen s215098
=======================================================================================

WHAT THE ABBREVIATIONS ARE
  NNBT       N,N-bis(2-hydroxypropyl)-p-toluidine. The laccase substrate. Laccase
             oxidises the C-H alpha to the tertiary nitrogen; the hemiaminal collapses
             and releases LACTALDEHYDE (2-hydroxypropanal) plus the secondary amine.
             That released lactaldehyde is what Purpald quantifies, which is why the
             standard curve is lactaldehyde. Dissolved in acetonitrile - see
             REACTION_MIX_FLOW_SCALE.
  Purpald    4-amino-3-hydrazino-5-mercapto-1,2,4-triazole in strong NaOH. Adds to an
             aldehyde, then the adduct must be OXIDISED BY DISSOLVED O2 to become the
             purple chromophore read at A530. Air-saturated buffer holds only ~250 uM
             O2, so aldehyde much above that gives the unoxidised intermediate -
             YELLOW-BROWN, not purple - and the response stops being monotonic.
  guaiacol   2-methoxyphenol, the redox mediator arm. 1 mM final in the reaction.
  lactaldehyde  the Purpald standard, tube D6. See LAC_STOCK_UM / LAC_GRADIENT_UM.

One dilution series feeds three readouts. Everything you change per run is a RUNTIME
PARAMETER in the Opentrons app - you should not need to open this file.

  NNBT plate  (slot 1, Heater-Shaker)  two arms: without and with guaiacol
  ABTS plate  (slot 5, swapped in)     activity, read immediately at A414/A734
  MALDI target(slot 5, swapped again)  spotted during the NNBT incubation
  readout plate (slot 5, swapped last) the 10x dilution of the NNBT plate, see below
  tube rack   (slot 7, MALDI only)     back on deck with the CAPPED matrix tube, while
                                       the 300 uL tips wait off-deck until Purpald

ABTS PLATE - there is NO ABTS-radical standard curve. The positive control is PaDa-1
at 1:1000, made up by hand into ONE tube (there is no neat PaDa-1 powder to run):

       1      2     3     4          col 1   PaDa-1 1:1000, a vertical triplicate,
  A  PaDa    NC1   NC1   NC1                 and the ONLY column that takes the
  B  1:1000  NC2   NC2   NC2                 PaDa mix
  C   |      Enz1  Enz1  Enz1        cols 2+ NC1, NC2, then the enzymes from row C,
  D   .      Enz2  ...                       horizontal triplicates, laccase mix
  E   .      Enz3                            (enzymes 7-8 spill into cols 5-7)
  F   .      ...up to Enz6
  G   .       .                      A dot is a well with reaction mix but no sample:
  H   .       .                      the mix goes in by whole columns (see below).

NNBT PLATE - 32 triplicate groups, 8 per 3-column block
  cols  1-3   NC2, NC3, 4 lactaldehyde standards, enz 9-10   mix: NNBT (NC3: no-NNBT)
  cols  4-6   the same 8, with guaiacol                      mix: guaiacol (NC3: no-NNBT+gua)
  cols  7-9   enzymes 1-8                                    mix: NNBT
  cols 10-12  enzymes 1-8, with guaiacol                     mix: guaiacol
  -> 10 enzymes maximum. NC2 = heat-inactivated enzyme,
     NC3 = active enzyme in reaction mix WITHOUT NNBT.

NC1 IS NOT ON THIS PLATE ANY MORE. A control block is 8 rows and 10 enzymes need two
of them, so two control slots had to go: NC1, and the fifth lactaldehyde point. NC1
was the one to drop because the 2026-09-15 run measured NC1 = 0.948 and NC2 = 1.008,
a difference of one plate-noise SD - NC2 is NC1 plus protein, so it is the better
blank and it subsumes NC1. NC3 stayed because it is the only well that gives the
reagent+plate optical floor (0.610), which every chromophore number is measured from.
NC1 SURVIVES ON THE ABTS PLATE, out of a dilution well in column 3 (see build_layout).

WELL RECIPE   10 uL spike + 140 uL reaction mix = 150 uL, 2 h @ 40 C.
              Then the plate is read TWICE, see THE 10x READOUT below.

THE 10x READOUT   Purpald's purple needs the aldehyde to stay in a narrow window. The
2026-09-15 run put every lactaldehyde standard past it: 500-6000 uM read 0.94-1.12
against a 0.929 blank, r2 = 0.03 over a 12-fold range - a flat line, because the
response peaks somewhere below 500 uM and collapses to brown above it. So after the
incubation every well is diluted 10x into milliQ on a SECOND plate, and BOTH plates
are developed and read:
  undiluted plate  135 uL left + 45 uL Purpald = 180 uL   (25% reagent, as before)
  10x plate         15 uL sample + 135 uL milliQ + 50 uL Purpald = 200 uL
The 10x plate is the one whose standards should be on scale. The undiluted plate keeps
continuity with every previous run, and comparing the two says which side of the peak
a sample sat on: 10x ~ undiluted/10 means the linear limb, anything higher means it
was browning. Both develop ON the Heater-Shaker, 10 min at 40 C / 1000 rpm, because
the purple step is an O2 oxidation and the shaking is what gets O2 back in - so they
are done one after the other, not side by side.

MULTI-DISPENSE   wherever one reagent goes to several wells in a row, ONE aspirate
serves several of them, with a touch-tip after every dispense: spikes 2 wells per
p20 load, reaction mix 2 wells per p300 load, Purpald 5 columns per load. Each p300
load carries a 20 uL disposal volume that is blown into the trash - build_layout()
has already added that to the volumes it asks you to pour. TWO PLACES DO NOT DO THIS:
the ABTS mix, because 2 x 190 uL will not fit a p300, and everything on the MALDI
target, where a touch-tip is forbidden (see the spotting section).

TWO WAYS TO SPIKE   the NNBT plate is spiked 8-CHANNEL: a dilution column holds eight
samples in the same row order as the destination block, so one aspirate loads all eight
nozzles. The ABTS plate cannot do that - its PaDa triplicate runs DOWN a column and its
enzymes start at row C - so every ABTS spike is SINGLE-NOZZLE, one sample at a time.
That is 4+n tips instead of two whole columns, so it is cheaper, and it is slower only
before the mix goes in, which is where slowness does not matter.

THE ABTS MIX STILL GOES IN BY WHOLE COLUMNS, 8-channel, because ABTS is kinetic and 36
single-nozzle dispenses would smear the start times across the plate. A column is all
one mix or nothing, which is why every PaDa well lives in column 1. The price is that
the unused wells of a part-filled column get reaction mix with no sample in them. They
are blank, nobody reads them, and the mix is cheap.

SLOT 5 IS USED THREE TIMES   tube rack -> ABTS plate -> MALDI target. They are never
needed at once. It is the middle deck column, so BOTH pipettes can reach every well of
everything that lands there - see rule 3 below.

RUN ORDER
  1  buffer into the dilution wells
  2  enzyme / NC3 / heat-inactivated stocks in, mixed      -> all at one molar conc
  3  lactaldehyde standard curve built by serial dilution
  4  PaDa-1 1:1000 moved from its tube into a dilution well      (rack is about to go)
  5  10 uL spikes into the NNBT plate, all four blocks           - 8-channel
  6  PAUSE - swap the tube rack for the ABTS plate
  7  10 uL spikes into the ABTS plate                            - SINGLE-NOZZLE
  8  reaction mixes into the NNBT plate
       blocks 1-2 one well at a time  (NC3 needs a different mix from its neighbours)
       blocks 3-4 8-channel by column (all one mix)
  9  ABTS mixes - LAST, because ABTS is kinetic and starts on contact.
       column 1 gets the PaDa mix (H2O2), every other column the laccase mix
 10  PAUSE - ABTS plate to the reader, seal the NNBT plate, MALDI target in slot 5,
       300 uL tips out of slot 7, tube rack + capped matrix tube into slot 7
 11  incubate 2 h @ 40 C, pausing at each MALDI timepoint (MALDI_TIMEPOINTS_MIN,
       t = 0 / 30 / 120 min) to unseal, premix + spot, reseal. ONE row of the target
       per sample and time running across its columns; the matrix tube is opened and
       closed ONCE PER SESSION, not per spot
 12  PAUSE - MALDI target out, EMPTY READOUT PLATE into slot 5, 300 uL tips back into
       slot 7, and a FRESH 20 uL rack into slot 3 (the 10x pass uses a whole rack)
 13  135 uL milliQ into the readout plate, then 15 uL out of every NNBT well into it,
       8-channel, fresh tips per column, mixed
 14  45 uL Purpald into the UNDILUTED plate on the shaker; develop 10 min. It is
       quenched first on purpose: it still holds active enzyme at full NNBT, while the
       10x plate is already running ten times slower
 15  PAUSE - undiluted plate to the reader (save as "-undiluted"), 10x plate onto the
       Heater-Shaker
 16  shake the 10x plate 2 min at 1000 rpm BEFORE Purpald - 15 uL under 135 uL will
       not mix itself (diffusion moves ~0.9 mm in the 12 min it sat; the column is
       ~4 mm), and a concentrated bolus meeting Purpald browns and stays brown
 17  50 uL Purpald into the 10x plate, develop, read A530 (save as "-10x")

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

MAX_ENZYMES = 10                              # 8 in the enzyme block + 2 that spill
                                              # into the control block's spare rows.
                                              # check_block_budget() keeps this in
                                              # step with LAC_GRADIENT_UM.


def add_parameters(p):
    """Everything you set per run. Shown in the app at run setup."""
    p.add_int(display_name='Number of enzymes', variable_name='n_enz', default=3,
              minimum=1, maximum=MAX_ENZYMES,
              description=f'Max {MAX_ENZYMES}: each enzyme needs a well in two blocks.')
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
    p.add_bool(display_name='MALDI spotting', variable_name='maldi_on', default=False,
               description='Spot the target during the NNBT incubation.')
    p.add_str(display_name='MALDI start row', variable_name='maldi_row', default='A',
              choices=[{'display_name': r, 'value': r} for r in 'ABCDEFGHIJKLMNOP'],
              description='First target row. ONE ROW PER SAMPLE: NC2, then the '
                          'enzymes. Timepoints run across the columns.')


metadata = {
    'apiLevel': '2.20',
    'protocolName': 'W-O-SC-10 - NNBT/guaiacol + ABTS + MALDI',
    'description': (
        'Dilutes every enzyme to one molar concentration, plates an NNBT/Purpald assay '
        'with and without guaiacol plus three negative controls and a lactaldehyde '
        'curve, spikes an ABTS activity plate from the same dilutions, and spots a '
        'MALDI target at set timepoints during the incubation.'),
}


# =======================================================================================
# SECTION 1 - THINGS THAT RARELY CHANGE   (everything per-run is a parameter above)
# =======================================================================================
ENZYME_NAMES = ['Lac-01', 'Lac-02', 'Lac-03', 'Lac-04', 'Lac-05',
                'Lac-06', 'Lac-07', 'Lac-08', 'Lac-09', 'Lac-10']   # plate-map labels
ENZ_TARGET_AUTO_FRACTION = 0.90        # auto target = this x the weakest enzyme, in uM

# Lactaldehyde positive control. FIVE points - the sixth was dropped to make room for
# NC3, so the controls fill one 3-column block exactly.
LAC_STOCK_UM = 100_000.0              # tube D6. 1 M diluted 1:10 = 100000 uM
LAC_GRADIENT_UM = [400, 1200, 2500, 5000]         # uM IN THE 150 uL REACTION
# FOUR points, not five: the fifth slot went to the 10th enzyme (see the header).
# The numbers come off the 2026-09-15 run, where 500-6000 uM was one flat line on the
# blank. With the 10x readout a well at C uM in the reaction is read at C x 15/200 =
# C/13.3, so these four are read at 30, 90, 188 and 375 uM - spanning the window the
# old run showed the response lives in, with the top point near the collapse so the
# curve shows you where its own ceiling is. Stock draws from the 100 mM D6 tube are
# 12, 36, 75 and 150 uL, all well clear of the p20 minimum, so the tube is unchanged.
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

# ABTS POSITIVE CONTROL - PaDa-1 at 1:1000, made up BY HAND into one tube. There is no
# ABTS-radical standard curve: we have no such standard, so the positive control is a
# known-active enzyme instead of a concentration series. There is no neat PaDa-1 powder,
# so only the 1:1000 strength runs. To bring a second strength back, add it to both
# lists - the ABTS column and dilution column 3 size themselves from PADA_NAMES.
PADA_NAMES = ['PaDa-1 1:1000']
PADA_TUBES = ['A5']
PADA_FILL_UL = 100.0                   # PRE-MADE and transferred, not diluted, so these
                                       # wells do NOT hold DILUTION_WELL_VOL_UL. 30 uL
                                       # is spiked out in two p20 loads (20 then 10) and
                                       # the tip must stay under the surface for both.
                                       # 40 left the tip 1.35 mm ABOVE the liquid. 60
                                       # cleared the last draw by 0.01 mm, which is not
                                       # clearance. 80 gave 0.56 mm, and that is already
                                       # in the regime where ASPIRATE_MIN_HEIGHT_MM has
                                       # overridden the half-depth rule. 100 leaves
                                       # 70 uL and 0.95 mm at the worst point, and 20 uL
                                       # more of a hand-made solution costs nothing.
                                       # The NC wells used to need a short fill like this
                                       # too, for the 8-channel column aspirate. They no
                                       # longer do - see build_layout().

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
PURPALD_VOL_UL = 50.0                  # 10x READOUT PLATE: 150 + 50 -> 200 uL
PURPALD_UNDIL_VOL_UL = 45.0            # UNDILUTED PLATE: the 10x pass has taken 15 uL
                                       # out, so it holds 135 uL, not 150. 45 keeps the
                                       # reagent fraction at the same 25% - same NaOH,
                                       # same chemistry - at the cost of 180 uL instead
                                       # of 200, i.e. ~10% less pathlength than older
                                       # runs. A constant factor across every well, so
                                       # nothing WITHIN the plate is affected. Using
                                       # 50 here instead would push it to 27% reagent.
ABTS_MIX_VOL_UL = 190.0                # 10 + 190 = 200 uL on the ABTS plate
REACTION_VOL_UL = SPIKE_VOL_UL + REACTION_MIX_VOL_UL

# --- the 10x readout -----------------------------------------------------------------
# 15 uL of a finished 150 uL reaction into 135 uL of milliQ, on a second plate, then
# Purpald as usual. milliQ and not buffer on purpose: the Purpald reagent's NaOH
# outweighs the assay buffer 10-20x either way, so it sets the development pH
# regardless, and water simply gives it less to fight.
READOUT_SAMPLE_UL = 15.0
READOUT_WATER_UL = 135.0
READOUT_VOL_UL = READOUT_SAMPLE_UL + READOUT_WATER_UL      # 150, as in the reaction
READOUT_FACTOR = READOUT_VOL_UL / READOUT_SAMPLE_UL        # 10x
READOUT_DRAW_MM = 2.0                  # into the finished NNBT well. The MALDI-spotted
                                       # wells are down to ~135 uL by now (5 uL a round
                                       # out of replicate 1 only), which is still ~4 mm
                                       # of liquid over a 2 mm tip.
READOUT_DISPENSE_MM = 1.0              # UNDER the 135 uL of water, not onto it
READOUT_MIX_REPS = 3
READOUT_MIX_UL = 20.0                  # the p20's whole barrel
READOUT_MIX_RPM = 1000                 # the shake that actually does the mixing, run
READOUT_MIX_MIN = 2                    # BEFORE Purpald - see the header. The p20
                                       # strokes above only break the layer up; 20 uL
                                       # in a 150 uL well cannot homogenise it.

REPLICATES = 3                         # triplicates, side by side in a row
WELLS_PER_COLUMN = 8
ROW_LETTERS = list('ABCDEFGH')
BLOCK_COLS = REPLICATES                # a 3-column block holds 8 triplicate groups
CTRL_FIXED_SLOTS = 2                   # NC2 and NC3, the two controls that stayed
CTRL_SPILL = max(0, MAX_ENZYMES - WELLS_PER_COLUMN)   # enzymes 9-10, in the control
                                       # block's last rows. check_block_budget().

INCUBATION_TEMP_C = 40                 # inside the module's 37-95 C range
INCUBATION_RPM = 250
INCUBATION_MIN = 120 #??? should be 120
DEVELOP_TEMP_C = 40
DEVELOP_RPM = 1000                     # not 3000: parafilm off, wells brim-full
DEVELOP_MIN = 10
REACTION_MIX_FLOW_SCALE = 0.5          # NNBT is in acetonitrile and drips at full speed
# Pipetting speed. The defaults (7.6 / 94 uL/s) make mixing crawl; these are your
# SC-08 values. Blow-out stays at each pipette's own default - a fast blow-out into a
# shallow well spatters and foams protein.
P20_FLOW_UL_S = 20.0                   # default 7.6, hardware max 24
P300_FLOW_UL_S = 150.0                 # default 94, hardware max 275
P20_BLOWOUT_UL_S = 10                 # p20 default: gentle
P300_BLOWOUT_UL_S = 110.0               # p300 default: gentle for this pipette
MALDI_FLOW_UL_S = 7.6                  # spotting only. A 1 uL droplet placed onto a
                                       # flat steel target at 15 uL/s splashes.

# MALDI: 5 uL sample + 20 uL milliQ = 1:5, one spot per condition. The 1:5 dilutions are
# made on the SAME dilution plate, from column MALDI_DIL_FIRST_COL on - columns 1-3 hold
# the assay dilutions, so there is no second plate to swap in.
#
# PREMIX. Matrix is NOT spotted on its own any more: each premix well gets
#   20 uL milliQ + 25 uL matrix + 5 uL sample = 50 uL, and 2 uL of that is spotted.
# That is the same 1:1 sample:matrix spot the old 1 uL matrix + 1 uL sample made on
# the target, just mixed in a well. It fixes three failures seen on the robot:
#   - 1 uL of acetonitrile matrix (the p20 minimum) was often not picked up or not
#     dispensed. Acetonitrile vapour in a dry tip pushes a 1 uL plug straight out.
#     25 uL, a pre-wetted tip, a slow aspirate and a pause fix that.
#   - 25 uL on a 6.9 mm FLAT floor is ~0.7 mm deep and creeps to the wall, leaving the
#     centre - where the tip is - dry. 50 uL with acetonitrile in it wets the floor.
#   - the 2 uL spot draw itself came up EMPTY: the tip mixed, arced up out of the well
#     and back down, and aspirated air. See MALDI_PREMIX_TIP_MM and the spotting loop -
#     the tip now goes deeper and never leaves the liquid between the mix and the draw.
MALDI_SAMPLE_UL = 5.0
MALDI_WATER_UL = 20.0
MALDI_MATRIX_UL = 25.0                 # from a capped tube, see MATRIX_TUBES;
                                       # drawn as 20 + 5, the tip holds 20
MALDI_PREMIX_UL = MALDI_WATER_UL + MALDI_MATRIX_UL + MALDI_SAMPLE_UL
MALDI_MATRIX_PREWET_REPS = 2           # saturate the tip with acetonitrile vapour first
MALDI_MATRIX_DELAY_S = 1.0             # after the matrix aspirate, before moving
MALDI_MATRIX_DISPENSE_MM = 3.0         # above the water in the premix well: the matrix
                                       # tip never touches liquid it could carry back
# HOW DEEP THE TIP GOES IN THE PREMIX WELL. 50 uL on a 6.86 mm flat floor is only
# 1.35 mm deep, so this number is almost all of the submersion there is. It used to be
# 0.8 mm, which left 0.34 mm under the tip at the bottom of a mix stroke - inside the
# stack of labware-offset error (+-0.3 mm is normal), the concave meniscus an
# acetonitrile mix pulls in the CENTRE of the well (where the tip is), and the
# acetonitrile that boils off while the well waits its turn. The tip came back up dry
# and spotted nothing. 0.5 mm is still 5x the p20 orifice bore, so it cannot seal on
# the floor, and it doubles the margin. check_maldi_geometry() enforces both ends.
MALDI_PREMIX_TIP_MM = 0.5              # mix and draw height in the 50 uL premix well
MALDI_PREMIX_MIN_TIP_MM = 0.3          # below this a p20 tip starts to seal on a flat
                                       # floor and moves nothing
# Water and matrix go into EVERY premix well of a round before the first sample does,
# so the last well of a round sits open for ~15 min with 25 uL of acetonitrile in it.
# It loses roughly 0.5-1 uL/min. The geometry check below must pass for a well that has
# already lost this much, not for a full one.
MALDI_PREMIX_EVAP_UL = 10.0
MALDI_PREMIX_MIX_REPS = 5
MALDI_PREMIX_MIX_UL = 8.0              # a bigger stroke drops the surface below the tip
                                       # - check_maldi_geometry() enforces it
MALDI_SPOT_UL = 2.0                    # premix per spot
MALDI_SPOT_SURPLUS_UL = 1.0            # drawn extra, never dispensed: no air on steel
MALDI_SPOT_DRAW_DELAY_S = 1.0          # hold the tip still under the surface after the
                                       # spot draw: 3 uL of a half-acetonitrile mix does
                                       # not finish entering the tip the instant the
                                       # plunger stops, and lifting early leaves it in
                                       # the well
MALDI_MIX_REPS = 3                     # aspirate/dispense on the target spreads the drop
MALDI_MIX_UL = 1.0                     # the p20's minimum, half the 2 uL spot
MALDI_SPOT_HEIGHT_MM = 0.2             # from the WELL BOTTOM, and the well is only
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

# WHEN TO SPOT. Not an interval any more - the useful timepoints are not evenly spaced.
# t=0 is the baseline, t=30 catches the fast phase, and t=120 is the endpoint; the
# intermediate points the old even-interval scheme spent matrix on were already
# collected on an earlier run and added nothing.
MALDI_TIMEPOINTS_MIN = [0, 30, 120]

# TARGET LAYOUT: ONE ROW PER SAMPLE, TIME RUNNING ACROSS.
#   rows     NC2 heat-inactivated, then the enzymes in enzyme-number order
#   columns  t0 plain, t0 guaiacol, t30 plain, t30 guaiacol, t120 plain, t120 guaiacol
# So a row is one sample's whole time course with both arms side by side, which is how
# you actually want to read it off the instrument. The old layout put one ROUND on two
# rows, which split a sample's time course across the target and - once enzymes 9-10
# moved into the control block - listed the samples out of order.
MALDI_ARMS = 2                         # plain, then guaiacol, within each timepoint


def maldi_spot(row_index, round_index, arm_index, start_row):
    """(row letter, 1-based column) of one spot."""
    return (MALDI_ROWS[start_row + row_index],
            MALDI_ARMS * round_index + arm_index + 1)
# What gets spotted: the heat-inactivated control (NC2) and the enzymes, one spot each
# per arm. NC1 (buffer) and NC3 (no NNBT) add nothing to the mass spectrum.
def maldi_per_arm(n_enz):
    return 1 + n_enz                   # NC2 + the enzymes: one TARGET ROW each


def maldi_spot_count(n_enz):
    """Total spots = rows x timepoints x arms. Also the number of premix wells."""
    return maldi_per_arm(n_enz) * len(MALDI_TIMEPOINTS_MIN) * MALDI_ARMS


def maldi_rows(layout):
    """Target row labels, top to bottom: NC2 then the enzymes in enzyme-number
    order - which is NOT plate-block order, since enzymes 9-10 sit in the control
    block."""
    return ([layout['heat_dil']['name']]
            + [d['name'] for d in layout['enz_dils']])


def matrix_plan(spots):
    """[(tube, uL of matrix the run draws from it)], in the order the tubes are used.
    Pour each tube with MATRIX_TUBE_DEAD_UL on top."""
    need = spots * MALDI_MATRIX_UL
    per_tube = MATRIX_TUBE_MAX_UL - MATRIX_TUBE_DEAD_UL
    per_tube -= per_tube % MALDI_MATRIX_UL                # whole draws only
    plan = []
    for tube in MATRIX_TUBES:
        if need <= 1e-9:
            break
        plan.append((tube, min(need, per_tube)))
        need -= plan[-1][1]
    if need > 1e-9:
        raise ValueError(f'MALDI needs {spots * MALDI_MATRIX_UL:.0f} uL of matrix but '
                         f'{len(MATRIX_TUBES)} tube(s) hold {len(MATRIX_TUBES) * per_tube:.0f}'
                         ' uL. Add a tube to MATRIX_TUBES or drop a timepoint from '
                         'MALDI_TIMEPOINTS_MIN.')
    return plan


MALDI_DIL_FIRST_COL = 4                # first dilution-plate column for the 1:5 wells
MALDI_DIL_WELLS = (12 - MALDI_DIL_FIRST_COL + 1) * 8   # columns 4-12 = 72 wells


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
SLOT_DILUTION = '6'                    # dilution plate; MALDI 1:5 wells from col 4
SLOT_TIPRACK_300 = '7'                 # 300 uL tips                 (north: 10, empty)
SLOT_MALDI_TUBERACK = '7'              # MALDI ONLY: the tube rack comes back here with
                                       # the capped matrix tube while the 300 uL tips
                                       # (idle until Purpald) wait off-deck. Nothing
                                       # else works: 5 holds the target, 2 is east of
                                       # the Heater-Shaker (max 53 mm, the rack is 54),
                                       # 4 and 9 are north of single-nozzle targets,
                                       # and 10/11 are out of reach for the H1 nozzle.
SLOT_RESERVOIR = '8'                   # reagents                    (north: 11, empty)
# Slots 2, 4, 9, 10, 11 MUST STAY EMPTY - clearance, not spare space.

NNBT_PLATE = 'eppendorf_96_wellplate_350ul'
READOUT_PLATE = NNBT_PLATE             # THE SAME PLATE TYPE AS THE NNBT PLATE, and not
                                       # by preference: the 10x plate finishes the run
                                       # ON the Heater-Shaker, and
                                       # opentrons_96_flat_bottom_adapter refuses to
                                       # stack corning_96_wellplate_360ul_flat.
                                       # Matching plates also means matching well
                                       # geometry, so the two A530 reads differ only by
                                       # the dilution and the Purpald volume.
DILUTION_PLATE = 'corning_96_wellplate_360ul_flat'
ABTS_PLATE = 'corning_96_wellplate_360ul_flat'
MALDI_PLATE = 'maldi_384_wellplate'    # Lukas's real definition, face at 18.0 mm
RESERVOIR_LOADNAME = 'nest_12_reservoir_15ml'
TUBERACK = '3d_printed_tuberack_1.5ml'
HS_ADAPTER = 'opentrons_96_flat_bottom_adapter'  # MUST be declared or every Z is wrong

TUBERACK_ROWS = 4                      # rack is 4 rows (A-D) x 6 columns = 24 tubes
TUBE_HEAT_INACT = 'C6'                 # fixed positions, so they never move
TUBE_LACTALDEHYDE = 'D6'
MATRIX_TUBES = ['B5', 'B6']            # MALDI matrix (DHB/acetonitrile), CAPPED, in the
                                       # rack's free columns 5-6. The second tube is only
                                       # used when the first cannot cover the run.
MATRIX_TUBE_MAX_UL = 1500.0
MATRIX_TUBE_DEAD_UL = 150.0            # the conical tip of the tube the p20 cannot empty

# Reservoir. Four different reaction mixes, all premixed off-deck by you.
RES_NNBT = ['A1', 'A2']                # buffer + NNBT                (2 wells: >13 mL)
RES_BUFFER = 'A3'                      # dilution buffer + the NC1 spike
RES_PURPALD = ['A4']
RES_ABTS = ['A5']                      # laccase ABTS mix. One well: even at 10 enzymes
                                       # the plate takes 6 columns x 8 x 190 = 9.1 mL,
                                       # and A6 is now the PaDa mix.
RES_PADA = ['A6']                      # PaDa-1 ABTS mix: H2O2 instead of Cu2+, and a
                                       # different pH. PaDa-1 is a peroxygenase, so it
                                       # cannot run in the laccase mix. Column 1 of the
                                       # ABTS plate is the only column it touches.
RES_GUAIACOL = ['A7', 'A8']            # guaiacol buffer + NNBT + guaiacol
RES_NO_NNBT = ['A9']                   # NC3: reaction mix WITHOUT NNBT
RES_NO_NNBT_GUA = ['A10']              # NC3 guaiacol: no NNBT, with guaiacol
RES_WATER = ['A11', 'A12']             # milliQ: the MALDI 1:5 dilutions AND the 10x
                                       # readout. TWO WELLS, and that is not optional -
                                       # the readout alone is 135 uL x 96 = 13.0 mL and
                                       # a well only carries RES_USABLE_UL once its
                                       # RES_DEAD_UL is in. spread() sizes it; both
                                       # MALDI and the readout draw through water_src()
                                       # in run(), which empties A11 before A12.
# A12 used to be free. The MALDI matrix used to live there, but acetonitrile evaporates
# out of an open reservoir well over a 2 h run - it is in a capped tube now
# (MATRIX_TUBES), which is what freed A12 for this.
RES_USABLE_UL = 13_000.0               # 15 mL nominal, minus fill margin
RES_DEAD_UL = 1_000.0                  # pour this much extra so tips never hit air

P20_MIN, P20_MAX = 1.0, 20.0
P300_MIN, P300_MAX = 20.0, 200.0

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

# The ABTS plate's two mixes. A COLUMN TAKES ONE OR THE OTHER, NEVER BOTH - the mix is
# dispensed 8-channel by column to keep the kinetic start times together, and an
# 8-channel dispense cannot tell rows apart. Column 1 is therefore the PaDa column.
MIX_ABTS, MIX_PADA = 'abts', 'pada'
ABTS_MIX_RESERVOIR = {MIX_ABTS: RES_ABTS, MIX_PADA: RES_PADA}
ABTS_MIX_LABEL = {MIX_ABTS: 'laccase ABTS mix', MIX_PADA: 'PaDa-1 ABTS mix (H2O2)'}
PADA_COLUMN = 1                        # the only ABTS column that takes the PaDa mix


def _tube(index):
    """Tube-rack well for the index-th enzyme (0-based): 4 rows per column."""
    return f'{ROW_LETTERS[index % TUBERACK_ROWS]}{index // TUBERACK_ROWS + 1}'


def _dil(index):
    """Dilution-plate well for the index-th dilution (0-based), column-then-row."""
    return f'{ROW_LETTERS[index % WELLS_PER_COLUMN]}{index // WELLS_PER_COLUMN + 1}'


def _abts_row_slots():
    """(row, first column) of each horizontal triplicate on the ABTS plate, in fill
    order. Column 1 is the PaDa column, so the rows start at column 2; each block is
    REPLICATES columns wide and the blocks run left to right without end. The first two
    slots are the negative controls, which is what puts the enzymes on row C."""
    col = PADA_COLUMN + 1
    while True:
        for row in range(WELLS_PER_COLUMN):
            yield row, col
        col += REPLICATES


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

    # ---- dilution plate -------------------------------------------------------------
    #   col 1  NC2, NC3, lac x4, enz 9-10 -> NNBT control blocks, both arms
    #   col 2  enzymes 1-8                -> NNBT enzyme blocks, and the ABTS enzymes
    #   col 3  PaDa-1 1:1000, then NC1    -> the ABTS plate ONLY
    #   col 4+ empty until the incubation -> the MALDI 1:5 dilutions
    #
    # COLUMNS 1 AND 2 ARE ROW-ALIGNED WITH THE NNBT BLOCKS. Every NNBT block is 8 groups,
    # one per ROW, REPLICATES columns wide, so a column of a block is "one replicate of
    # all 8 samples" - exactly what an 8-channel dispenses. One aspirate is ONE Z for all
    # eight nozzles, so every well in such a column must hold the SAME volume and be
    # drawn down on the same schedule.
    #
    # COLUMN 3 IS NOT, and does not have to be: the ABTS plate is spiked single-nozzle,
    # so its two wells are free to carry PADA_FILL_UL instead of a full well.
    #
    # NC1 LIVES IN COLUMN 3 NOW. It is off the NNBT plate entirely (see the header), so
    # its only job is the ABTS plate's buffer control - 30 uL, single-nozzle. Column 3 is
    # the one column no 8-channel pass ever aspirates, so it does not have to be row-
    # aligned with anything and its wells are free to carry PADA_FILL_UL instead of a
    # full 200. That is the same 30-out-of-100 draw the PaDa well has already proved.
    #
    # NC2 is still drawn single-nozzle out of its 200 uL col-1 well for ABTS: 30 uL out
    # of 200 leaves 1.49 mm of clearance on the worst draw. It is safe ONLY because no
    # 8-channel pass touches column 1 after the ABTS spikes have pulled that well below
    # its neighbours - see step 7 in run().
    dilutions = []
    warnings = []                                     # samples that went in neat
    V_MAIN, V_PADA = DILUTION_WELL_VOL_UL, PADA_FILL_UL

    def add_dil(well, name, kind, address, stock_conc, fill, note='', conc_to=None):
        """stock_conc None = a buffer-only well (NC1, and the fillers)."""
        conc = None
        if stock_conc is None:
            stock, buf = 0.0, fill
        elif conc_to is None and stock_conc < target:
            # TOO DILUTE FOR THE TARGET: go in NEAT instead of refusing to build. Only
            # samples take this path (conc_to None) - a lactaldehyde standard that
            # missed its concentration would corrupt the curve, so those still raise.
            stock, buf, conc = fill, 0.0, stock_conc
            warnings.append(
                f'{name} (tube {address}) is only {stock_conc:.3g} uM, below the '
                f'{target:.3g} uM target - dispensed NEAT at {stock_conc:.3g} uM '
                f'({stock_conc / target:.0%} of target). Normalise its rates by hand.')
            note = (note + '; ' if note else '') + f'NEAT - {stock_conc:.3g} uM, not target'
        else:
            stock, buf = dilution_recipe(stock_conc, conc_to or target, total=fill)
            conc = conc_to or target
        d = {'name': name, 'well': well, 'src_kind': kind, 'src': address,
             'stock_conc': stock_conc or 0.0, 'stock_vol': stock, 'buffer_vol': buf,
             'fill': fill, 'note': note, 'conc': conc}
        dilutions.append(d)
        return d

    def _w(row, col):
        return f'{ROW_LETTERS[row]}{col}'

    def filler(well):
        """Idle nozzles would otherwise aspirate air. Buffer costs nothing and the
        NNBT rows these feed are unused anyway, so the spikes land in wells nobody
        reads."""
        return add_dil(well, f'(unused row {well[0]})', 'res', RES_BUFFER, None, V_MAIN,
                       'filler: keeps the 8th nozzle out of a dry well')

    # --- column 1: the NNBT control block, in block row order ------------------------
    heat_dil = add_dil(_w(0, 1), 'NC2 heat-inactivated', 'tube', TUBE_HEAT_INACT,
                       heat_um, V_MAIN)
    nc3_dil = add_dil(_w(1, 1), f'NC3 ({enz_um[nc3_index]["name"]})', 'tube',
                      _tube(nc3_index), enz_um[nc3_index]['conc'], V_MAIN,
                      'second aliquot, keeps NC3 off the enzyme well so neither is '
                      'over-drawn')
    spike_factor = REACTION_VOL_UL / SPIKE_VOL_UL          # 15x
    lac_dils, src_conc, src = [], float(LAC_STOCK_UM), ('tube', TUBE_LACTALDEHYDE)
    for row, um in enumerate(sorted(LAC_GRADIENT_UM, reverse=True),
                             start=CTRL_FIXED_SLOTS):
        want = um * spike_factor
        d = add_dil(_w(row, 1), f'Lac-std {um:g} uM', src[0], src[1], src_conc, V_MAIN,
                    'serial' if src[0] == 'dil' else '', conc_to=want)
        lac_dils.append(d)
        if LAC_SERIAL:                                 # chain: next point comes from
            src_conc, src = want, ('dil', d['well'])   # this one. See LAC_SERIAL.

    # The control block's last CTRL_SPILL rows: the enzymes that do not fit the 8-row
    # enzyme block. They are ordinary enzyme samples in an ordinary NNBT/guaiacol mix -
    # the only thing "control block" means for them is which column they sit in.
    spill_dils, spill_enz = [], []
    for k in range(CTRL_SPILL):
        row, i = CTRL_FIXED_SLOTS + len(LAC_GRADIENT_UM) + k, WELLS_PER_COLUMN + k
        if i < n:
            d = add_dil(_w(row, 1), enz_um[i]['name'], 'tube', _tube(i),
                        enz_um[i]['conc'], V_MAIN,
                        'enzyme 9-10: the control block is where they fit')
            spill_enz.append(d)
        else:
            d = filler(_w(row, 1))
        spill_dils.append(d)

    # --- column 2: enzymes 1-8, in enzyme block row order ----------------------------
    head = min(n, WELLS_PER_COLUMN)
    block_dils = [add_dil(_w(i, 2), e['name'], 'tube', _tube(i), e['conc'], V_MAIN)
                  for i, e in enumerate(enz_um[:head])]
    fillers = [filler(_w(i, 2)) for i in range(head, WELLS_PER_COLUMN)]
    enz_dils = block_dils + spill_enz        # every enzyme, in enzyme-number order

    # --- column 3: the ABTS-only wells -----------------------------------------------
    # PaDa-1 is made up by hand as powder in buffer and transferred whole - no dilution
    # on deck. NC1 follows it: buffer only, and only the ABTS plate reads it.
    pada_dils = []
    for row, (name, tube) in enumerate(zip(PADA_NAMES, PADA_TUBES)):
        d = {'name': name, 'well': _w(row, 3), 'src_kind': 'tube', 'src': tube,
             'stock_conc': 0.0, 'stock_vol': 0.0, 'buffer_vol': 0.0,
             'fill': V_PADA, 'note': 'pre-made by hand, transferred not diluted'}
        dilutions.append(d)
        pada_dils.append(d)
    nc1_dil = add_dil(_w(len(PADA_NAMES), 3), 'NC1 buffer only', 'res', RES_BUFFER,
                      None, V_PADA,
                      'ABTS ONLY - NC1 is off the NNBT plate. Column 3 is never '
                      'aspirated 8-channel, so it needs no row alignment')

    # ---- NNBT plate: four blocks of eight groups ------------------------------------
    def ctrl_block(block, arm, mix, nomix):
        """One control block: NC2, NC3, the lactaldehyde standards, then whichever
        enzymes spilled out of the enzyme block. Row order matches dilution column 1,
        which is what lets one 8-channel aspirate load all eight."""
        groups = [
            {'label': heat_dil['name'], 'arm': arm, 'mix': mix,
             'src': ('dil', heat_dil['well'])},
            {'label': nc3_dil['name'], 'arm': arm, 'mix': nomix,
             'src': ('dil', nc3_dil['well'])},
        ] + [{'label': d['name'], 'arm': arm, 'mix': mix, 'src': ('dil', d['well'])}
             for d in lac_dils + spill_dils]
        for slot, g in enumerate(groups):
            g['wells'] = _block_wells(block, slot)
        return groups

    def enz_block(block, arm, mix):
        groups = [{'label': d['name'], 'arm': arm, 'mix': mix,
                   'src': ('dil', d['well'])} for d in block_dils + fillers]
        for slot, g in enumerate(groups):
            g['wells'] = _block_wells(block, slot)
        return groups

    nnbt = (ctrl_block(0, ARM_PLAIN, MIX_NNBT, MIX_NO_NNBT)
            + ctrl_block(1, ARM_GUA, MIX_GUA, MIX_NO_NNBT_GUA)
            + enz_block(2, ARM_PLAIN, MIX_NNBT)
            + enz_block(3, ARM_GUA, MIX_GUA))

    # ---- ABTS plate -----------------------------------------------------------------
    # Column 1: the PaDa-1 strength(s) as VERTICAL triplicates, so the whole positive
    # control sits in the one column that takes the H2O2 mix. Everything else is a
    # HORIZONTAL triplicate from column 2 on: NC1, NC2, then the enzymes from row C.
    # The fillers are not here - they exist only to keep an 8-channel column aspirate off
    # a dry well, and no 8-channel pass reads this plate's sources.
    if len(PADA_NAMES) * REPLICATES > WELLS_PER_COLUMN:
        raise ValueError(f'{REPLICATES} replicates x {len(PADA_NAMES)} PaDa strengths '
                         f'needs {len(PADA_NAMES) * REPLICATES} rows, but a column has '
                         f'{WELLS_PER_COLUMN}.')
    if max(int(d['well'][1:]) for d in dilutions) >= MALDI_DIL_FIRST_COL:
        raise ValueError('the assay dilutions reach into the MALDI dilution columns - '
                         'raise MALDI_DIL_FIRST_COL.')
    abts = []
    for i, d in enumerate(pada_dils):
        row0 = i * REPLICATES
        abts.append({'label': d['name'], 'src': ('dil', d['well']), 'mix': MIX_PADA,
                     'wells': [f'{ROW_LETTERS[row0 + k]}{PADA_COLUMN}'
                               for k in range(REPLICATES)]})
    slots = _abts_row_slots()
    for d in [nc1_dil, heat_dil] + enz_dils:
        row, col = next(slots)
        abts.append({'label': d['name'], 'src': ('dil', d['well']), 'mix': MIX_ABTS,
                     'wells': [f'{ROW_LETTERS[row]}{col + k}' for k in range(REPLICATES)]})
    abts_columns = max(int(w[1:]) for g in abts for w in g['wells'])
    # Which mix each whole column takes. Only column 1 is PaDa; the rest are laccase.
    abts_col_mix = {c: (MIX_PADA if c == PADA_COLUMN else MIX_ABTS)
                    for c in range(1, abts_columns + 1)}

    # ---- reagent volumes ------------------------------------------------------------
    # Blocks 1-2 are filled one well at a time (NC3 needs its own mix), so they cost
    # only the wells actually used. Blocks 3-4 are 8-channel, so a part-filled column
    # still costs a full column.
    per_mix = {}
    ctrl_wells = {}
    for g in nnbt[:16]:                                   # the two control blocks
        per_mix[g['mix']] = per_mix.get(g['mix'], 0.0) + REACTION_MIX_VOL_UL * REPLICATES
        ctrl_wells[g['mix']] = ctrl_wells.get(g['mix'], 0) + REPLICATES
    enz_cols = BLOCK_COLS if head else 0
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

    nnbt_columns = 2 * BLOCK_COLS + (2 * BLOCK_COLS if head else 0)

    # PURPALD GOES ONTO TWO PLATES: PURPALD_UNDIL_VOL_UL into the NNBT plate once the
    # 10x pass has drawn 15 uL out of it, and PURPALD_VOL_UL into the readout plate.
    # Both are 8-channel by column, so each costs a whole column of wells and each
    # load's disposal volume is paid by all eight nozzles.
    purpald_total = sum(
        vol * WELLS_PER_COLUMN * nnbt_columns
        + P300_MIN * p300_loads(nnbt_columns, vol) * WELLS_PER_COLUMN
        for vol in (PURPALD_VOL_UL, PURPALD_UNDIL_VOL_UL))

    # milliQ: the 10x readout's diluent, plus the MALDI premix water if it is on. The
    # readout dispense is one column per load, so a disposal volume would buy nothing
    # and multi_dispense() is called with disposal=0 - hence no waste term.
    #
    # THE SLACK TERM IS NOT PADDING. The readout draws milliQ one whole 8-channel load
    # at a time - READOUT_WATER_UL x 8 = 1080 uL - and an aspirate cannot be split
    # across two reservoir wells. So each well can only serve floor(its share / 1080)
    # columns and strands the remainder, and spread()'s even split has no idea. Without
    # one spare load per well the 10-enzyme + MALDI run dies eleven columns in with
    # 'milliQ ran out: 1080 uL wanted, 660 uL left' - found by opentrons_simulate,
    # which is the only thing that would have found it.
    readout_load = READOUT_WATER_UL * WELLS_PER_COLUMN
    readout_water_total = (readout_load * nnbt_columns
                           + readout_load * len(RES_WATER))

    buffer_total = (sum(d['buffer_vol'] for d in dilutions)
                    + SPIKE_VOL_UL * REPLICATES)         # NC1 on the ABTS plate
    # The ABTS mixes go in by WHOLE COLUMNS, so a part-filled column still costs eight
    # wells of mix. That waste is the price of keeping the kinetic starts together.
    abts_per_mix = {}
    for key in abts_col_mix.values():
        abts_per_mix[key] = (abts_per_mix.get(key, 0.0)
                             + ABTS_MIX_VOL_UL * WELLS_PER_COLUMN)
    return {
        'target': target, 'enzymes': enz_um, 'dilutions': dilutions,
        'enz_dils': enz_dils, 'nc3_dil': nc3_dil, 'heat_dil': heat_dil,
        'lac_dils': lac_dils, 'pada_dils': pada_dils, 'nc1_dil': nc1_dil,
        'nnbt': nnbt, 'abts': abts,
        'readout_water_total': readout_water_total,
        'nnbt_columns': nnbt_columns,
        'abts_columns': abts_columns, 'abts_col_mix': abts_col_mix,
        'per_mix': per_mix, 'buffer_total': buffer_total,
        'purpald_total': purpald_total,
        'abts_per_mix': abts_per_mix,
        'warnings': warnings,
    }


def p300_loads(n_dests, vol, disposal=P300_MIN):
    """How many p300 aspirates multi_dispense() needs for n_dests wells of `vol`."""
    if not n_dests:
        return 0
    per_load = max(1, int((P300_MAX - disposal) // vol))
    return math.ceil(n_dests / per_load)


def spread(total_ul, wells, reagent):
    """Split a reagent evenly over the fewest reservoir wells that hold it.

    EVERY well carries its own RES_DEAD_UL on top of its share, so the test is
    share + dead <= usable, not total <= usable. Splitting therefore does not halve
    the requirement - it adds a dead volume. Getting this wrong is how the 10x
    readout's 13.0 mL of milliQ came out as "one well" and then asked for 14.0 mL to
    be poured into a well that holds 13.0."""
    need = 1
    while (total_ul / need + RES_DEAD_UL > RES_USABLE_UL) and need < len(wells):
        need += 1
    if total_ul / need + RES_DEAD_UL > RES_USABLE_UL:
        raise ValueError(
            f'{reagent} needs {total_ul / 1000:.1f} mL plus {RES_DEAD_UL / 1000:g} mL '
            f'dead volume per well, which will not fit {len(wells)} well(s) of '
            f'{RES_USABLE_UL / 1000:g} mL. Configure another reservoir well.')
    return {w: total_ul / need for w in wells[:need]}


# --- preflight: refuse to load on a rule we have already been burned by ----------------
LABWARE_HEIGHT_MM = {'opentrons_96_tiprack_20ul': 64.7,
                     'opentrons_96_tiprack_300ul': 64.5, TUBERACK: 54.0,
                     RESERVOIR_LOADNAME: 31.4, NNBT_PLATE: 14.3, MALDI_PLATE: 18.0,
                     'TRASH': 999.0}


def check_deck():
    """Slot 9 can never be a single-nozzle target (trash in 12 is north of it) and
    nothing tall may sit north of one. opentrons_simulate does NOT catch the trash
    case - it passed the layout the robot then refused. Checked for the deck as set up
    AND for the deck during MALDI spotting, when slots 5 and 7 hold something else."""
    setup = {str(HS_SLOT): NNBT_PLATE, SLOT_SWAP: TUBERACK,
             SLOT_TIPRACK_20: 'opentrons_96_tiprack_20ul',
             SLOT_DILUTION: DILUTION_PLATE,
             SLOT_TIPRACK_300: 'opentrons_96_tiprack_300ul',
             SLOT_RESERVOIR: RESERVOIR_LOADNAME, '12': 'TRASH'}
    maldi = dict(setup)
    maldi[SLOT_SWAP] = MALDI_PLATE
    del maldi[SLOT_TIPRACK_300]                   # tips wait off-deck ...
    maldi[SLOT_MALDI_TUBERACK] = TUBERACK         # ... while the matrix rack is on
    readout = dict(setup)                         # the 10x pass: tube rack gone, the
    readout[SLOT_SWAP] = READOUT_PLATE            # readout plate in its place
    problems = []
    for phase, occupied in (('setup', setup), ('MALDI', maldi), ('readout', readout)):
        for slot, what in occupied.items():
            if slot == '12':
                continue
            north = str(int(slot) + 3)
            if north == '12':
                problems.append(f'[{phase}] {what} in slot {slot} has the fixed trash '
                                'bin north of it - this is the 2026-09-04 failure.')
            elif north in occupied and LABWARE_HEIGHT_MM.get(occupied[north], 0) >= 50:
                problems.append(f'[{phase}] {what} in slot {slot} has {occupied[north]} '
                                'north of it - too tall for a single-nozzle move.')
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


def check_maldi_geometry():
    """The premix well is shallow. The p20 has to stay under the surface at the bottom
    of every mix stroke and after the spot draw - that is where it pooled at the wall
    and came up dry.

    Checked against the WORST well of a round, not a fresh one: MALDI_PREMIX_EVAP_UL of
    acetonitrile is already gone by the time the last well is sampled."""
    if MALDI_PREMIX_TIP_MM < MALDI_PREMIX_MIN_TIP_MM:
        raise ValueError(f'MALDI premix tip at {MALDI_PREMIX_TIP_MM:g} mm above a FLAT '
                         f'floor - a p20 tip seals below ~{MALDI_PREMIX_MIN_TIP_MM:g} mm '
                         'and moves nothing.')
    area = math.pi * (6.86 / 2) ** 2
    full = MALDI_PREMIX_UL - MALDI_PREMIX_EVAP_UL
    for what, left in (('premix mix stroke', full - MALDI_PREMIX_MIX_UL),
                       ('spot draw', full - MALDI_SPOT_UL - MALDI_SPOT_SURPLUS_UL)):
        surface = left / area
        if MALDI_PREMIX_TIP_MM > surface - ASPIRATE_MARGIN_MM:
            raise ValueError(
                f'MALDI {what}: a well that has lost {MALDI_PREMIX_EVAP_UL:g} uL to '
                f'evaporation holds {left:.1f} uL, surface at {surface:.2f} mm, tip at '
                f'{MALDI_PREMIX_TIP_MM:g} mm - it would draw AIR. Lower '
                'MALDI_PREMIX_TIP_MM, shrink MALDI_PREMIX_MIX_UL, or raise the premix '
                'volume.')


check_deck()
def check_block_budget():
    """A control block is one column of 8 rows and it has to hold NC2, NC3, every
    lactaldehyde standard AND whatever enzymes do not fit the 8-row enzyme block.
    Change LAC_GRADIENT_UM or MAX_ENZYMES without changing the other and the plate
    silently stops adding up - so it raises at import instead."""
    spill = max(0, MAX_ENZYMES - WELLS_PER_COLUMN)
    used = CTRL_FIXED_SLOTS + len(LAC_GRADIENT_UM) + spill
    if used != WELLS_PER_COLUMN:
        raise ValueError(
            f'a control block has {WELLS_PER_COLUMN} rows but the layout wants {used}: '
            f'{CTRL_FIXED_SLOTS} fixed controls (NC2, NC3) + {len(LAC_GRADIENT_UM)} '
            f'lactaldehyde standards + {spill} enzyme(s) spilling out of the enzyme '
            f'block. Add or drop a standard, or change MAX_ENZYMES.')


def check_maldi_timepoints():
    """The timepoints have to be a rising list inside the incubation, and the two
    columns each one takes have to exist on the target."""
    t = MALDI_TIMEPOINTS_MIN
    if not t or list(t) != sorted(set(t)):
        raise ValueError(f'MALDI_TIMEPOINTS_MIN must be a rising list of distinct '
                         f'minutes, got {t}.')
    if t[0] < 0 or t[-1] > INCUBATION_MIN:
        raise ValueError(f'MALDI_TIMEPOINTS_MIN {t} runs outside the '
                         f'{INCUBATION_MIN} min incubation.')
    if MALDI_ARMS * len(t) > MALDI_COLS:
        raise ValueError(f'{len(t)} timepoints x {MALDI_ARMS} arms needs '
                         f'{MALDI_ARMS * len(t)} target columns; there are {MALDI_COLS}.')


check_mix_geometry()
check_maldi_geometry()
check_block_budget()
check_maldi_timepoints()


# Colours for the Opentrons app's labware map. The app shows a coloured, named liquid
# in every tube, reservoir well and plate well, so the setup screen doubles as the
# "what goes where" sheet.
COLOUR = {
    'enzyme':    '#2f7d5f',            # green  - the enzymes under test
    'nc1':       '#78706f',            # grey   - buffer only
    'nc2':       '#3f6d8f',            # blue   - heat-inactivated
    'nc3':       '#7b4fa0',            # purple - active enzyme, no NNBT
    'lac':       '#c2621a',            # amber  - lactaldehyde standards
    'pada':      '#1b7f8f',            # teal   - PaDa-1 positive control
    'nnbt':      '#8e24aa',            # violet - NNBT reaction mix
    'gua':       '#b8860b',            # gold   - guaiacol reaction mix
    'no_nnbt':   '#c2185b',            # pink   - reaction mix without NNBT
    'buffer':    '#0288d1',            # cyan   - assay buffer
    'purpald':   '#c62828',            # red    - Purpald
    'abts_mix':  '#00695c',            # dark teal - laccase ABTS mix
    'pada_mix':  '#00838f',            # cyan-teal - PaDa-1 ABTS mix (H2O2)
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
    if name.startswith('PaDa'):
        return COLOUR['pada']
    return COLOUR['enzyme']


# =======================================================================================
# SECTION 5 - REPORTING   (run log and `python3 <this file>`)
# =======================================================================================
SYMBOL = {MIX_NNBT: 'N', MIX_GUA: 'G', MIX_NO_NNBT: 'x', MIX_NO_NNBT_GUA: 'y'}


def render(layout, n_enz, maldi_on=False, maldi_row='A'):
    """Everything you need to set the deck up and read the plates afterwards."""
    out = ['=' * 78,
           f'  {n_enz} enzymes x 2 arms (plain + guaiacol), NC2 + NC3 on the plate '
           f'(NC1 is ABTS-only), {len(LAC_GRADIENT_UM)} lactaldehyde standards',
           f'  READ TWICE: the NNBT plate itself, and its {READOUT_FACTOR:g}x dilution '
           'on the readout plate',
           '  ABTS positive control: ' + ' and '.join(PADA_NAMES) + ' (no ABTS standard)',
           f'  every enzyme diluted to {layout["target"]:.2f} uM',
           '=' * 78, '']
    if layout.get('warnings'):
        out += ['!' * 78, '  WARNING - NOT EVERY SAMPLE IS AT THE TARGET CONCENTRATION'] \
               + [f'  - {w}' for w in layout['warnings']] + ['!' * 78, '']
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
    # Irregular on purpose: the PaDa positive control runs DOWN column 1, everything
    # else runs ACROSS from column 2. A grid is the only readable way to show that.
    apos = {w: i + 1 for i, g in enumerate(layout['abts']) for w in g['wells']}
    cols = layout['abts_columns']
    out += ['', '-' * 78,
            f'  ABTS PLATE (slot {SLOT_SWAP}, swapped in after the dilutions)',
            '  P=PaDa mix (H2O2)  A=laccase ABTS mix  --=mix but NO sample, do not read',
            '       ' + ''.join(f'{c:>5}' for c in range(1, cols + 1))]
    for r in ROW_LETTERS:
        cells = []
        for c in range(1, cols + 1):
            sym = 'P' if layout['abts_col_mix'][c] == MIX_PADA else 'A'
            w = f'{r}{c}'
            cells.append((f'{sym}{apos[w]:02d}' if w in apos else f'{sym}--').rjust(5))
        out.append(f'  {r}  ' + ''.join(cells))
    out += ['', '   #  wells          content']
    for i, g in enumerate(layout['abts']):
        out.append(f'  {i + 1:02d}  {",".join(g["wells"]):<14} {g["label"]} '
                   f'[{ABTS_MIX_LABEL[g["mix"]]}]')

    # --- dilution plate ---
    out += ['', '-' * 78,
            f'  DILUTION PLATE (slot 6), {DILUTION_WELL_VOL_UL:g} uL per well',
            '  well  from        stock uL  buffer uL  content']
    for d in layout['dilutions']:
        origin = {'tube': 'tube ', 'res': 'res  '}.get(d['src_kind'], 'dil  ') + d['src']
        out.append(f'  {d["well"]:<5} {origin:<11} {d["stock_vol"]:>8.1f} '
                   f'{d["buffer_vol"]:>10.1f}  {d["name"]}'
                   f'{"  [" + d["note"] + "]" if d["note"] else ""}')

    # --- the 10x readout plate ---
    out += ['', '-' * 78,
            f'  READOUT PLATE (slot {SLOT_SWAP}, swapped in AFTER the incubation)',
            f'  Use a {READOUT_PLATE} - the SAME type as the NNBT plate. It finishes '
            'on the Heater-Shaker, whose adapter will not take a Corning flat plate.',
            f'  Same 96-well map as the NNBT plate, well for well: '
            f'{READOUT_WATER_UL:g} uL milliQ first, then {READOUT_SAMPLE_UL:g} uL out '
            f'of the matching NNBT well = {READOUT_FACTOR:g}x.',
            f'  Purpald: {PURPALD_UNDIL_VOL_UL:g} uL into the NNBT plate (it is down '
            f'to {REACTION_VOL_UL - READOUT_SAMPLE_UL:g} uL by then), '
            f'{PURPALD_VOL_UL:g} uL into this one.',
            '  Both develop on the Heater-Shaker, one after the other. Read A530 twice:',
            '    "-undiluted"  the NNBT plate   - continuity with earlier runs',
            f'    "-10x"        this plate      - the standards should be on scale here']

    # --- what to pour ---
    needs = {RES_BUFFER: ('assay buffer', layout['buffer_total'])}
    for key, total in layout['per_mix'].items():
        for w, v in spread(total, MIX_RESERVOIR[key], MIX_LABEL[key]).items():
            needs[w] = (MIX_LABEL[key], v)
    for w, v in spread(layout['purpald_total'], RES_PURPALD, 'Purpald').items():
        needs[w] = ('Purpald reagent (both plates)', v)
    for key, total in layout['abts_per_mix'].items():
        for w, v in spread(total, ABTS_MIX_RESERVOIR[key], ABTS_MIX_LABEL[key]).items():
            needs[w] = (ABTS_MIX_LABEL[key], v)
    water_total = layout['readout_water_total']
    if maldi_on:
        spots = maldi_spot_count(n_enz)
        water_total += MALDI_WATER_UL * spots
    for w, v in spread(water_total, RES_WATER, 'milliQ').items():
        needs[w] = ('milliQ (10x readout' + (' + MALDI)' if maldi_on else ')'), v)
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
        tubes[d['src']] = (vol + (d['stock_vol'] or d['fill']), names + [d['name']])
    for tube in sorted(tubes, key=lambda t: (int(t[1:]), t[0])):
        vol, names = tubes[tube]
        out.append(f'    {tube:<4} {" + ".join(names):<34} >= {vol + 300:.0f} uL')
    if maldi_on:
        out.append(f'  MATRIX - CAPPED tube(s), put in the rack when it goes back into '
                   f'slot {SLOT_MALDI_TUBERACK} for MALDI')
        for tube, vol in matrix_plan(spots):
            out.append(f'    {tube:<4} {"MALDI matrix (DHB/acetonitrile)":<34} '
                       f'{vol + MATRIX_TUBE_DEAD_UL:.0f} uL')

    # --- MALDI ---
    if maldi_on:
        rows = maldi_rows(layout)
        start = MALDI_ROWS.index(maldi_row)
        if start + len(rows) > len(MALDI_ROWS):
            raise ValueError(f'MALDI needs {len(rows)} rows from {maldi_row}; the '
                             f'target has {len(MALDI_ROWS)}. Start higher up.')
        out += ['', '-' * 78,
                f'  MALDI - ONE ROW PER SAMPLE, time running across the columns',
                f'  {len(MALDI_TIMEPOINTS_MIN)} timepoints at t = '
                + ', '.join(f'{t:g}' for t in MALDI_TIMEPOINTS_MIN) + ' min, '
                f'{maldi_spot_count(n_enz)} spots in total',
                f'  premix per spot: {MALDI_WATER_UL:g} uL milliQ + {MALDI_MATRIX_UL:g} '
                f'uL matrix + {MALDI_SAMPLE_UL:g} uL sample, {MALDI_SPOT_UL:g} uL spotted',
                '', '  target map (row = sample, column = timepoint x arm)', '']
        pad, cell = max(len(r) for r in rows) + 2, 7
        out.append(' ' * (2 + pad)
                   + ''.join(f't{t:g} min'.center(cell * MALDI_ARMS)
                             for t in MALDI_TIMEPOINTS_MIN))
        out.append(' ' * (2 + pad)
                   + ''.join(('plain' if a == 0 else 'gua').rjust(cell)
                             for _ in MALDI_TIMEPOINTS_MIN
                             for a in range(MALDI_ARMS)))
        for i, label in enumerate(rows):
            cells = ''.join(
                '{}{}'.format(*maldi_spot(i, r, a, start)).rjust(cell)
                for r in range(len(MALDI_TIMEPOINTS_MIN)) for a in range(MALDI_ARMS))
            out.append(f'  {label:<{pad}}' + cells)
        out.append('')
        out += ['  plain = NNBT mix, gua = guaiacol mix. Only NC2 and the enzymes are',
                '  spotted - NC1 and NC3 add nothing to the mass spectrum.',
                '',
                '  The matrix tube is opened and closed ONCE PER SPOTTING SESSION, not',
                '  once per spot: the matrix tip is dispensed above the water and never',
                '  touches it, and every premix well holds the same thing at that point,',
                '  so one tip serves the whole round.']
    return out


# =======================================================================================
# SECTION 6 - THE PROTOCOL
# =======================================================================================
def run(protocol):
    prm = protocol.params
    n = prm.n_enz
    rounds = len(MALDI_TIMEPOINTS_MIN) if prm.maldi_on else 0
    maldi_spots = maldi_spot_count(n) if prm.maldi_on else 0

    # ---- the batch, straight from the app ------------------------------------------
    batch = [{'name': ENZYME_NAMES[i], 'mg_ml': getattr(prm, f'mg_ml_{i + 1}'),
              'mw_kda': getattr(prm, f'mw_{i + 1}')} for i in range(n)]
    heat_um = to_um(prm.heat_mg_ml, prm.heat_mw)
    layout = build_layout(batch, heat_um, prm.target_um, prm.nc3_enz - 1)

    # ---- MALDI fits? check now, not two hours into the incubation -------------------
    # BEFORE render() and the liquid map, both of which call matrix_plan(), and 'add a
    # tube to MATRIX_TUBES' is the wrong advice when the real problem is that the run
    # does not fit the dilution plate or the target at all.
    if prm.maldi_on:
        if maldi_spots > MALDI_DIL_WELLS:
            raise ValueError(
                f'MALDI needs one premix well per spot = {maldi_spots}, but '
                f'dilution-plate columns {MALDI_DIL_FIRST_COL}-12 hold '
                f'{MALDI_DIL_WELLS}. Drop a timepoint from MALDI_TIMEPOINTS_MIN or run '
                'fewer enzymes.')
        if MALDI_ROWS.index(prm.maldi_row) + maldi_per_arm(n) > len(MALDI_ROWS):
            raise ValueError(
                f'MALDI needs {maldi_per_arm(n)} rows from {prm.maldi_row} (NC2 + '
                f'{n} enzymes, one row each); the target has {len(MALDI_ROWS)}. '
                'Start higher up.')
        matrix_plan(maldi_spots)                       # raises if the tubes are short

    for e in layout['enzymes']:                       # the conversion, on the record
        protocol.comment(f'  {e["name"]}: {e["mg_ml"]:g} mg/mL / {e["mw_kda"]:g} kDa '
                         f'= {e["conc"]:.2f} uM')
    for line in render(layout, n, prm.maldi_on, prm.maldi_row):
        protocol.comment(line)
    if layout['warnings']:                          # make the operator see it
        protocol.pause('WARNING: ' + ' | '.join(layout['warnings'])
                       + '  -  Resume to run anyway.')
    # ---- deck ----------------------------------------------------------------------
    hs = protocol.load_module('heaterShakerModuleV1', HS_SLOT)
    hs_adapter = hs.load_adapter(HS_ADAPTER)      # kept: the readout plate lands on it
    nnbt_plate = hs_adapter.load_labware(NNBT_PLATE)
    tuberack = protocol.load_labware(TUBERACK, SLOT_SWAP)
    dil_plate = protocol.load_labware(DILUTION_PLATE, SLOT_DILUTION)
    reservoir = protocol.load_labware(RESERVOIR_LOADNAME, SLOT_RESERVOIR)
    tr20 = protocol.load_labware('opentrons_96_tiprack_20ul', SLOT_TIPRACK_20)
    tr300 = protocol.load_labware('opentrons_96_filtertiprack_200ul', SLOT_TIPRACK_300)
    abts_plate = protocol.load_labware(ABTS_PLATE, protocol_api.OFF_DECK)
    maldi_target = (protocol.load_labware(MALDI_PLATE, protocol_api.OFF_DECK)
                    if prm.maldi_on else None)
    readout_plate = protocol.load_labware(READOUT_PLATE, protocol_api.OFF_DECK,
                                          label='10x readout plate')

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
            tube_needs[d['src']] = (vol + (d['stock_vol'] or d['fill']),
                                    names + [d['name']])
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
    for key, total in layout['abts_per_mix'].items():
        tag = 'pada_mix' if key == MIX_PADA else 'abts_mix'
        for w, v in spread(total, ABTS_MIX_RESERVOIR[key], ABTS_MIX_LABEL[key]).items():
            res_liquids.append((w, ABTS_MIX_LABEL[key].capitalize(), tag, v))
    water_total = layout['readout_water_total'] + MALDI_WATER_UL * maldi_spots
    water_plan = spread(water_total, RES_WATER, 'milliQ')
    for w, v in water_plan.items():
        res_liquids.append((w, 'milliQ (10x readout + MALDI)', 'water', v))
    if prm.maldi_on:
        spots = maldi_spots
        for tube, vol in matrix_plan(spots):
            tuberack[tube].load_liquid(liquid('MALDI matrix - CAPPED tube',
                                              COLOUR['matrix']), vol + MATRIX_TUBE_DEAD_UL)
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
            pip.reset_tipracks()                      # else the API still calls them used
        well = racks[key][nxt[key]]                   # row A when count == 8
        nxt[key] += count
        pip.pick_up_tip(well)

    # milliQ is the one reagent two different steps draw from, so the wells are
    # emptied in order through a single tracker rather than each step assuming a well.
    water_wells = sorted(water_plan, key=lambda w: int(w[1:]))
    water_left = dict(water_plan)

    def water_src(load):
        """A reservoir well with at least `load` uL of milliQ left in it. `load` is
        what the PIPETTE draws - 8 nozzles x the per-well volume for a column pass."""
        for w in water_wells:
            if water_left[w] >= load:
                water_left[w] -= load
                return reservoir[w]
        raise ValueError(
            f'milliQ ran out: {load:.0f} uL wanted, '
            f'{max(water_left.values()):.0f} uL left in the fullest of '
            f'{", ".join(water_wells)}. build_layout() and spread() are out of step.')

    def slow():
        """Gentle p20: droplets on the MALDI target, and acetonitrile matrix."""
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
    # Per round: one water tip, one matrix tip (it serves the whole round - see the
    # spotting loop), and one sample tip per spot.
    maldi_tips = rounds * (2 + MALDI_ARMS * maldi_per_arm(n))
    spike_cols = 2 if n else 1            # NNBT only: col 1, and col 2 if there are any
    abts_tips = len(layout['abts'])       # ABTS is single-nozzle: one tip per sample
    # The 10x readout is NOT in this number: it runs one fresh tip column per plate
    # column, which is a whole 96-tip rack on its own, and step 12 makes the operator
    # put a new rack in before it starts. Everything below has to fit the rack they
    # loaded at setup.
    tips20 = (len(layout['dilutions']) + spike_cols * WELLS_PER_COLUMN + abts_tips
              + maldi_tips + WELLS_PER_COLUMN)   # slack: pick20_column skips part-used
    have20 = len(racks['20']) - start['20']
    mid_run_racks = max(0, math.ceil((tips20 - have20) / len(racks['20'])))
    total_racks = 1 + mid_run_racks + 1            # loaded now, mid-run, then the 10x
    protocol.comment(f'  20 uL tips needed ~{tips20} before the 10x readout, available '
                     f'from {prm.tip20_row}{prm.tip20_col}: {have20}')
    protocol.comment(f'  *** HAVE {total_racks} FULL 20 uL RACKS READY: the one in slot '
                     f'{SLOT_TIPRACK_20} now, {mid_run_racks} more during the run, and '
                     'one at step 12 - the 10x readout spends a whole rack on its own, '
                     'one tip column per plate column. ***')
    if tips20 > have20:
        protocol.comment(f'  *** NOT ENOUGH 20 uL TIPS FOR THE RUN UP TO STEP 12: have '
                         f'{have20}, need ~{tips20}. It will PAUSE for a fresh rack '
                         f'{mid_run_racks} time(s) - have them open. ***')

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

    fill_of = {d['well']: d['fill'] for d in layout['dilutions']}

    def start_vol(well):
        """How much is in a dilution well before anything is spiked out of it."""
        return fill_of.get(well, DILUTION_WELL_VOL_UL)

    def pick20_column():
        while nxt['20'] % WELLS_PER_COLUMN:           # skip a part-used column
            nxt['20'] += 1
        take('20', p20, SLOT_TIPRACK_20, WELLS_PER_COLUMN)

    # HOW MUCH IS LEFT IN EVERY DILUTION WELL, tracked per well rather than per column.
    # It has to be per well now: the col-1 NC wells and the col-2 enzyme wells are drawn
    # 8-channel for the NNBT plate and then single-nozzle again for the ABTS plate, so
    # after the ABTS spikes the wells of a column no longer hold the same volume.
    # draw_at() reads this to place the tip, and raises rather than aspirate air.
    remaining = {d['well']: d['fill'] for d in layout['dilutions']}
    for d in layout['dilutions']:        # serial lac donations, if LAC_SERIAL is on,
        if d['src_kind'] == 'dil' and d['stock_vol'] > 0:    # left those wells short
            remaining[d['src']] -= d['stock_vol']

    def spike_column(src_col, dest_anchors, plate):
        """ONE 8-CHANNEL PASS. The dilution column holds eight different samples in
        the same row order as the destination block, so one aspirate loads all eight
        nozzles at once and each load is multi-dispensed across REPLICATES destination
        columns. All eight nozzles share one Z, which is only valid because every well
        in the column holds the same volume and is drawn down on the same schedule -
        so this may only run BEFORE the single-nozzle ABTS spikes, never after."""
        anchor = dil_plate[f'A{src_col}']
        col_wells = [w for w in (f'{r}{src_col}' for r in ROW_LETTERS) if w in remaining]

        def src(load):
            where = draw_at(anchor, remaining[f'A{src_col}'], load)
            for w in col_wells:                  # every nozzle draws the same load
                remaining[w] -= load
            return where

        pick20_column()
        multi_dispense(p20, SPIKE_VOL_UL, src, [plate[a] for a in dest_anchors],
                       height=SPIKE_HEIGHT_MM, touch_height=SPIKE_LOW_HEIGHT_MM,
                       disposal=0.0, touch_src=True)
        p20.drop_tip()

    def spike_single(src_well, dests, plate):
        """ONE SAMPLE, SINGLE-NOZZLE, into wells that need not share a row or a column.
        Still multi-dispensed: 2 wells per 20 uL p20 load, so a triplicate costs two
        aspirates (20 then 10). This is how the whole ABTS plate is spiked."""
        pick20()

        def src(load):
            where = draw_at(dil_plate[src_well], remaining[src_well], load)
            remaining[src_well] -= load
            return where

        multi_dispense(p20, SPIKE_VOL_UL, src, [plate[w] for w in dests],
                       height=SPIKE_HEIGHT_MM, touch_height=SPIKE_LOW_HEIGHT_MM,
                       disposal=0.0, touch_src=True)
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
    # 2/3/4  stocks in and mixed; lactaldehyde curve; PaDa-1 off the rack
    # ===================================================================================
    for d in layout['dilutions']:
        if d['stock_vol'] <= 0:                        # PaDa-1, handled below
            continue
        src = (tuberack[d['src']] if d['src_kind'] == 'tube'
               else draw_at(dil_plate[d['src']], DILUTION_WELL_VOL_UL, d['stock_vol']))
        transfer(d['stock_vol'], src, dil_well(d['well']),
                 mix_well=dil_plate[d['well']])

    # PaDa-1 is made up by hand and moved whole - no dilution, no mixing. It has to come
    # across now, while the tube rack is still on the deck: slot 5 is about to become the
    # ABTS plate.
    for d in layout['pada_dils']:
        transfer(d['fill'], tuberack[d['src']], dil_well(d['well']))

    # ===================================================================================
    # 5  spikes into the NNBT plate - all four blocks, 8-channel
    # ===================================================================================
    # THESE MUST COME FIRST. An 8-channel aspirate is one Z for all eight nozzles, which
    # only works while every well in the source column still holds the same volume. The
    # ABTS spikes below break that for A1/B1 and the enzyme wells, so nothing may
    # 8-channel out of columns 1 or 2 afterwards.
    p20.configure_nozzle_layout(style=ALL, tip_racks=[tr20])     # 8-channel for spikes
    ctrl_cols = [f'A{c}' for c in range(1, 2 * BLOCK_COLS + 1)]
    enz_cols_a = [f'A{c}' for c in range(2 * BLOCK_COLS + 1, 4 * BLOCK_COLS + 1)]
    spike_column(1, ctrl_cols, nnbt_plate)                   # controls, both arms
    if n:
        spike_column(2, enz_cols_a, nnbt_plate)              # enzymes, both arms
    p20.configure_nozzle_layout(style=SINGLE, start='H1', tip_racks=[tr20])

    # ===================================================================================
    # 6/7  tube rack out, ABTS plate in, spike it SINGLE-NOZZLE
    # ===================================================================================
    # One tip per sample, in layout order: PaDa-1 1:1000 down column 1, then
    # NC1, NC2 and the enzymes across from column 2. No 8-channel pass can do this - the
    # PaDa triplicates run down a column and the enzymes start at row C - and at 4+n tips
    # it is cheaper than the two tip columns the old aligned layout needed.
    protocol.move_labware(tuberack, protocol_api.OFF_DECK, use_gripper=False)
    protocol.move_labware(abts_plate, SLOT_SWAP, use_gripper=False)
    for g in layout['abts']:
        spike_single(g['src'][1], g['wells'], abts_plate)

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
    # 9  ABTS mixes - LAST liquid step: ABTS is kinetic and starts on contact
    # ===================================================================================
    # 8-CHANNEL, WHOLE COLUMNS, on purpose. 190 uL is more than half a p300, so there is
    # no multi-dispense to be had either way; going by column instead of by well is what
    # keeps 30-odd kinetic reactions from starting minutes apart. The cost is that the
    # empty rows of a part-filled column get mix with no sample - blank wells nobody
    # reads. That is also why the PaDa samples all live in column 1: a column takes ONE
    # mix, and an 8-channel dispense cannot tell rows apart.
    abts_res_of = {}
    for key, total in layout['abts_per_mix'].items():
        abts_res_of[key] = list(spread(total, ABTS_MIX_RESERVOIR[key],
                                       ABTS_MIX_LABEL[key]))[0]
    # ONE FRESH COLUMN OF TIPS PER MIX. A tip that has carried the PaDa mix carries H2O2,
    # and H2O2 in a laccase well is no longer a laccase assay. The laccase columns go
    # first, so the enzyme samples are the ones that wait the shortest.
    for key in (MIX_ABTS, MIX_PADA):
        cols = sorted(c for c, k in layout['abts_col_mix'].items() if k == key)
        if not cols:
            continue
        protocol.comment(f'  {ABTS_MIX_LABEL[key]} -> ABTS columns '
                         f'{", ".join(str(c) for c in cols)}')
        pick300_column()
        for c in cols:
            p300.aspirate(ABTS_MIX_VOL_UL, reservoir[abts_res_of[key]])
            p300.dispense(ABTS_MIX_VOL_UL,
                          abts_plate[f'A{c}'].bottom(z=REAGENT_HEIGHT_MM))
            p300.touch_tip()
        p300.drop_tip()

    # ===================================================================================
    # 10  hand over: ABTS plate to the reader, seal the NNBT plate, MALDI labware in
    # ===================================================================================
    protocol.pause(f'TAKE THE ABTS PLATE from slot {SLOT_SWAP} to the reader NOW '
                   '(A414/A734) - it is already reacting. Then resume.')
    protocol.move_labware(abts_plate, protocol_api.OFF_DECK, use_gripper=False)
    if prm.maldi_on:                       # the dilution plate STAYS: MALDI uses col 4+
        protocol.move_labware(maldi_target, SLOT_SWAP, use_gripper=False)
        # The matrix tube needs a rack on the deck. The 300 uL tips are idle until
        # Purpald, so they step out of slot 7 and the tube rack takes their place.
        protocol.pause(f'MALDI setup: take the 300 uL tip rack OUT of slot '
                       f'{SLOT_TIPRACK_300} and keep it clean - it goes back before '
                       'Purpald. Put the CAPPED matrix tube(s) '
                       f'{", ".join(t for t, _ in matrix_plan(maldi_spots))} in the tube '
                       f'rack and the rack in slot {SLOT_MALDI_TUBERACK}. Then resume.')
        protocol.move_labware(tr300, protocol_api.OFF_DECK, use_gripper=False)
        protocol.move_labware(tuberack, SLOT_MALDI_TUBERACK, use_gripper=False)
    protocol.pause(f'Seal the NNBT plate for the {INCUBATION_MIN} min incubation at '
                   f'{INCUBATION_TEMP_C} degC, then resume.')

    # ===================================================================================
    # 11  incubate, pausing at each MALDI timepoint to spot the target
    # ===================================================================================
    hs.set_and_wait_for_temperature(INCUBATION_TEMP_C)
    if not prm.maldi_on:
        hs.set_and_wait_for_shake_speed(INCUBATION_RPM)
        protocol.delay(minutes=INCUBATION_MIN)
        hs.deactivate_shaker()
    else:
        # One MALDI condition = one NNBT well, premixed with milliQ and matrix, one
        # spot. ONE TARGET ROW PER SAMPLE and time running across the columns - see
        # maldi_spot(). Only NC2 (the heat-inactivated control) and the enzymes are
        # spotted; NC1 and NC3 add nothing to the mass spectrum.
        start_row = MALDI_ROWS.index(prm.maldi_row)
        spotted = ([layout['heat_dil']['well']]
                   + [d['well'] for d in layout['enz_dils']])
        # Sorted by that list, NOT by plate order: enzymes 9-10 live in the control
        # block, so plate order would spot them between NC2 and enzyme 1. The spot key
        # is NC2, then the enzymes in enzyme-number order, on both arms.
        rank = {w: i for i, w in enumerate(spotted)}
        conds = [(g['label'], g['arm'], g['wells'][0], rank[g['src'][1]])
                 for g in layout['nnbt'] if g['src'][1] in spotted]  # replicate 1 only
        if len(conds) != 2 * maldi_per_arm(n):
            raise ValueError(f'MALDI picked {len(conds)} conditions, expected '
                             f'{2 * maldi_per_arm(n)} - maldi_per_arm() is out of step.')
        plain = sorted((c for c in conds if c[1] == ARM_PLAIN), key=lambda c: c[3])
        gua = sorted((c for c in conds if c[1] == ARM_GUA), key=lambda c: c[3])
        dil_i = [0]                     # running index into dilution-plate col 4 on
        maldi_offset = (MALDI_DIL_FIRST_COL - 1) * WELLS_PER_COLUMN

        def next_maldi_well():
            if dil_i[0] >= MALDI_DIL_WELLS:
                raise ValueError('the MALDI premix wells are full: '
                                 f'{maldi_spots} spots > {MALDI_DIL_WELLS} wells in '
                                 f'dilution columns {MALDI_DIL_FIRST_COL}-12. Drop a '
                                 'timepoint or run fewer enzymes.')
            w = _dil(maldi_offset + dil_i[0])
            dil_i[0] += 1
            return w

        # Which tube every matrix draw comes from, and what is left in each.
        plan = matrix_plan(maldi_spots)
        matrix_draws = [t for t, v in plan for _ in range(round(v / MALDI_MATRIX_UL))]
        matrix_left = {t: v + MATRIX_TUBE_DEAD_UL for t, v in plan}

        elapsed = 0.0                              # minutes of incubation done
        for r, t_min in enumerate(MALDI_TIMEPOINTS_MIN):
            if t_min > elapsed:                    # shake up to this timepoint
                hs.set_and_wait_for_shake_speed(INCUBATION_RPM)
                protocol.delay(minutes=t_min - elapsed)
                hs.deactivate_shaker()                 # never pipette a moving plate
                elapsed = t_min
            protocol.pause(f'MALDI round {r + 1}/{rounds} (t = {t_min:g} min): '
                           'UNSEAL the NNBT plate, then resume.')
            # A ROW PER SAMPLE, so the row index is the SAMPLE and the column is this
            # timepoint's pair of arms.
            todo = []                                  # (spot, NNBT well, premix well)
            for arm_i, arm in enumerate((plain, gua)):
                for i, (_, _, src_well, _) in enumerate(arm):
                    row, col = maldi_spot(i, r, arm_i, start_row)
                    todo.append((maldi_target[f'{row}{col}'], src_well,
                                 dil_plate[next_maldi_well()]))

            # WATER into every premix well of the round, ONE tip: the wells are empty,
            # so this tip only ever touches water and may go back to the reservoir.
            pick20()
            for _, _, well in todo:
                p20.aspirate(MALDI_WATER_UL, water_src(MALDI_WATER_UL))
                p20.dispense(MALDI_WATER_UL, well.bottom(z=DILUTION_HEIGHT_MM))
                clear(p20, well)
            p20.drop_tip()

            # MATRIX into every premix well of the round, ONE TIP AND ONE OPEN/CLOSE
            # PER SPOTTING SESSION - not per spot. It is safe for the same reason the
            # water pass is: the matrix is dispensed MALDI_MATRIX_DISPENSE_MM ABOVE the
            # water and never touches it, and at this moment every premix well holds
            # exactly the same thing - water and matrix, no sample yet, which only goes
            # in below on a clean tip. So there is nothing to cross-contaminate.
            # NO touch_tip on a premix well here, unlike everywhere else: 25 uL creeps
            # up the wall of a flat 6.9 mm well, and this tip goes straight back into
            # the matrix tube. A wet wall would carry water into the acetonitrile and
            # dilute the stock over the round's draws. Blow out above the liquid only.
            pick20()
            slow()                                     # acetonitrile drips at speed
            open_tube, prewet = None, True
            for _, _, well in todo:
                tube = matrix_draws.pop(0)
                if tube != open_tube:                  # only when the tube CHANGES
                    if open_tube is not None:
                        protocol.pause(f'CLOSE matrix tube {open_tube}, then resume.')
                    protocol.pause(
                        f'OPEN matrix tube {tube} (slot {SLOT_MALDI_TUBERACK}) for the '
                        f'whole of MALDI round {r + 1}, then resume.')
                    open_tube, prewet = tube, True
                # 25 uL does not fit a 20 uL tip: 20 + 5, both while the tube is open.
                done = 0.0
                while done < MALDI_MATRIX_UL - 1e-6:
                    stroke = min(MALDI_MATRIX_UL - done, P20_MAX)
                    src = draw_at(tuberack[tube], matrix_left[tube], stroke)
                    if prewet:              # saturate the tip with acetonitrile vapour
                        p20.mix(MALDI_MATRIX_PREWET_REPS, stroke, src)
                        prewet = False      # once per tube now, not once per spot
                    p20.aspirate(stroke, src)
                    matrix_left[tube] -= stroke
                    protocol.delay(seconds=MALDI_MATRIX_DELAY_S)
                    touch(p20, tuberack[tube])         # drop falls back into the tube
                    p20.dispense(stroke, well.bottom(z=MALDI_MATRIX_DISPENSE_MM))
                    done += stroke
                    if BLOW_OUT:
                        z = height(well, MALDI_WATER_UL + done) + BLOWOUT_ABOVE_MM
                        p20.blow_out(well.bottom(z=min(z, well.depth - 1.0)))
            fast()
            p20.drop_tip()
            protocol.pause(f'CLOSE matrix tube {open_tube}, then resume.')

            for spot, src_well, well in todo:
                # SAMPLE: a clean tip into the NNBT well, premix, spot 2 uL, spread it.
                pick20()
                premix = well.bottom(z=MALDI_PREMIX_TIP_MM)
                p20.aspirate(MALDI_SAMPLE_UL,
                             nnbt_plate[src_well].bottom(z=MALDI_DRAW_HEIGHT_MM))
                p20.dispense(MALDI_SAMPLE_UL, premix)
                # NO LOCATION on the mix or on the spot draw, deliberately. Passing one
                # makes the API arc the tip up to the top of the well and back down
                # between every command, and each of those exits drags premix out on the
                # tip, disturbs the shallow pool and gives the re-entry another chance to
                # land above the liquid - the tip came back down and aspirated nothing.
                # Without a location both act where the tip already is: under the
                # surface, from the dispense right through to the draw.
                p20.mix(MALDI_PREMIX_MIX_REPS, MALDI_PREMIX_MIX_UL)
                slow()                                 # from here the tip is on steel
                p20.aspirate(MALDI_SPOT_UL + MALDI_SPOT_SURPLUS_UL)
                protocol.delay(seconds=MALDI_SPOT_DRAW_DELAY_S)   # let it finish
                # NO blow-out, NO touch-tip, NO air gap on the target: air pushed
                # through the droplet sprays it onto neighbouring spots and leaves a
                # bubble, which is a hole in the crystal layer.
                p20.dispense(MALDI_SPOT_UL, spot.bottom(z=MALDI_SPOT_HEIGHT_MM))
                p20.mix(MALDI_MIX_REPS, MALDI_MIX_UL,     # 3 strokes spread the drop
                        spot.bottom(z=MALDI_SPOT_HEIGHT_MM))
                fast()
                p20.drop_tip()
            protocol.pause('RESEAL the NNBT plate, then resume.')
        if INCUBATION_MIN > elapsed:               # tail of the incubation, if any
            hs.set_and_wait_for_shake_speed(INCUBATION_RPM)
            protocol.delay(minutes=INCUBATION_MIN - elapsed)
            hs.deactivate_shaker()
        protocol.pause(f'Take the MALDI target from slot {SLOT_SWAP} and let the spots '
                       'dry - the matrix is already in them. Spot key is in this log. '
                       f'Take the tube rack out of slot {SLOT_MALDI_TUBERACK} and put the '
                       '300 uL tip rack back. Resume to finish.')
        protocol.move_labware(maldi_target, protocol_api.OFF_DECK, use_gripper=False)
        protocol.move_labware(tuberack, protocol_api.OFF_DECK, use_gripper=False)
        protocol.move_labware(tr300, SLOT_TIPRACK_300, use_gripper=False)

    # ===================================================================================
    # 12  the 10x readout plate goes on, and the 20 uL rack is replaced
    # ===================================================================================
    cols = [f'A{c}' for c in range(1, layout['nnbt_columns'] + 1)]
    protocol.pause(
        'Remove the seal from the NNBT plate. Put an EMPTY 96-well plate in slot '
        f'{SLOT_SWAP} (the 10x readout plate) and a FRESH 20 uL tip rack in slot '
        f'{SLOT_TIPRACK_20} - the 10x pass uses a whole rack, one column per plate '
        'column. Then resume.')
    protocol.move_labware(readout_plate, SLOT_SWAP, use_gripper=False)
    nxt['20'] = 0                                 # the operator just put a full rack in
    p20.reset_tipracks()                          # ... and the API has to be told

    # ===================================================================================
    # 13  135 uL milliQ, then 15 uL of every reaction on top of it
    # ===================================================================================
    # WATER FIRST, SAMPLE INTO IT. Dispensing 15 uL into standing water mixes far better
    # than dropping water onto 15 uL, and it rinses the tip on the way out.
    # 8-channel, one column per aspirate: 135 uL x 8 = 1080 uL leaves the reservoir per
    # load, which is what water_src() is charged. disposal=0 on purpose - one
    # destination per load means a disposal volume would be pure waste, and the full
    # blow-out it replaces is the more accurate ending anyway.
    pick300_column()
    for a in cols:
        p300.aspirate(READOUT_WATER_UL,
                      water_src(READOUT_WATER_UL * WELLS_PER_COLUMN))
        p300.dispense(READOUT_WATER_UL,
                      readout_plate[a].bottom(z=DILUTION_HEIGHT_MM))
        clear(p300, readout_plate[a], over=READOUT_WATER_UL)
    p300.drop_tip()

    # A FRESH TIP COLUMN PER PLATE COLUMN. These are 96 different samples; one carried
    # tip would cross-contaminate a whole row of the readout.
    p20.configure_nozzle_layout(style=ALL, tip_racks=[tr20])
    for a in cols:
        pick20_column()
        p20.aspirate(READOUT_SAMPLE_UL, nnbt_plate[a].bottom(z=READOUT_DRAW_MM))
        touch(p20, nnbt_plate[a])                 # drop falls back into the reaction
        p20.dispense(READOUT_SAMPLE_UL,
                     readout_plate[a].bottom(z=READOUT_DISPENSE_MM))
        # These strokes do NOT homogenise the well - 20 uL cannot stir 150. They break
        # the 15 uL up so it is not a discrete layer; the shake in step 16 is what
        # actually mixes it.
        p20.mix(READOUT_MIX_REPS, READOUT_MIX_UL,
                mix_at(readout_plate[a], READOUT_VOL_UL, READOUT_MIX_UL))
        clear(p20, readout_plate[a], over=READOUT_VOL_UL)
        p20.drop_tip()
    p20.configure_nozzle_layout(style=SINGLE, start='H1', tip_racks=[tr20])

    # ===================================================================================
    # 14  Purpald into the UNDILUTED plate, and develop it - this one first
    # ===================================================================================
    # ORDER MATTERS. The NNBT plate still holds active enzyme at full NNBT concentration
    # and is still making lactaldehyde; the readout plate is already running ten times
    # slower. Purpald's NaOH is the quench, so it goes to the plate that needs quenching.
    pur = list(spread(layout['purpald_total'], RES_PURPALD, 'Purpald'))[0]
    pick300_column()
    multi_dispense(p300, PURPALD_UNDIL_VOL_UL, reservoir[pur],   # several cols per load
                   [nnbt_plate[a] for a in cols])
    p300.drop_tip()

    hs.set_and_wait_for_temperature(DEVELOP_TEMP_C)
    hs.set_and_wait_for_shake_speed(DEVELOP_RPM)
    protocol.delay(minutes=DEVELOP_MIN)
    hs.deactivate_shaker()

    # ===================================================================================
    # 15  swap the plates on the Heater-Shaker
    # ===================================================================================
    # Both plates develop ON the shaker, one after the other, because the purple step is
    # an O2 oxidation and 1000 rpm is how O2 gets back into a 200 uL well. Developing one
    # of them static in slot 5 would not be the same assay.
    hs.open_labware_latch()
    protocol.pause(
        'READ THE UNDILUTED PLATE NOW: take it off the Heater-Shaker to the reader '
        '(A530) and save it as "-undiluted". Then put the 10x readout plate from slot '
        f'{SLOT_SWAP} onto the Heater-Shaker and resume.')
    protocol.move_labware(nnbt_plate, protocol_api.OFF_DECK, use_gripper=False)
    protocol.move_labware(readout_plate, hs_adapter, use_gripper=False)
    hs.close_labware_latch()

    # ===================================================================================
    # 16/17  mix the 10x plate, Purpald, develop
    # ===================================================================================
    # THE MIX THAT COUNTS. 15 uL under 135 uL does not mix itself: diffusion moves a
    # small molecule about 0.9 mm in the ~12 min this plate has been standing, and the
    # liquid column is ~4 mm. If Purpald met a concentrated bolus it would brown it, and
    # brown does not come back - which is exactly what happened to the standards above
    # 500 uM on 2026-09-15.
    hs.set_and_wait_for_temperature(DEVELOP_TEMP_C)
    hs.set_and_wait_for_shake_speed(READOUT_MIX_RPM)
    protocol.delay(minutes=READOUT_MIX_MIN)
    hs.deactivate_shaker()                        # never pipette a moving plate

    pick300_column()
    multi_dispense(p300, PURPALD_VOL_UL, reservoir[pur],
                   [readout_plate[a] for a in cols])
    p300.drop_tip()

    hs.set_and_wait_for_shake_speed(DEVELOP_RPM)
    protocol.delay(minutes=DEVELOP_MIN)
    hs.deactivate_shaker()
    hs.deactivate_heater()
    hs.open_labware_latch()
    protocol.comment(f'Done - read the 10x plate at A530 and save it as "-10x". The '
                     'undiluted plate was read in step 15. Plate maps, which are the '
                     'same for both, are at the top of this log.')


# =======================================================================================
# SECTION 7 - LOCAL PREVIEW   `python3 <this file>` - no robot, no opentrons needed
# =======================================================================================
if __name__ == '__main__':
    demo = [{'name': 'Lac-01', 'mg_ml': 2.40, 'mw_kda': 65.0},
            {'name': 'Lac-02', 'mg_ml': 1.85, 'mw_kda': 65.0},
            {'name': 'Lac-03', 'mg_ml': 3.10, 'mw_kda': 70.0}]
    print('\n'.join(render(build_layout(demo, to_um(2.40, 65.0)),
                           len(demo), True, 'A')))
