"""Builds the W-O-SC-10 operator setup guide PDF from the protocol's own logic.
Run: python3 build_wosc10_guide.py
"""
import importlib.util
import pathlib

_spec = importlib.util.spec_from_file_location(
    'wosc10', pathlib.Path(__file__).with_name(
        'Lac_W-O-SC-10_nnbt_guaiacol_abts_maldi_ot2_fixed (2).py'))
wosc10 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(wosc10)

build_layout, render, to_um = wosc10.build_layout, wosc10.render, wosc10.to_um
COLOUR = wosc10.COLOUR
PADA_NAMES, PADA_TUBES = wosc10.PADA_NAMES, wosc10.PADA_TUBES
TUBE_HEAT_INACT, TUBE_LACTALDEHYDE = wosc10.TUBE_HEAT_INACT, wosc10.TUBE_LACTALDEHYDE
RES_NNBT, RES_BUFFER, RES_PURPALD = wosc10.RES_NNBT, wosc10.RES_BUFFER, wosc10.RES_PURPALD
RES_ABTS, RES_GUAIACOL = wosc10.RES_ABTS, wosc10.RES_GUAIACOL
RES_PADA = wosc10.RES_PADA
RES_NO_NNBT, RES_NO_NNBT_GUA = wosc10.RES_NO_NNBT, wosc10.RES_NO_NNBT_GUA
RES_WATER, MATRIX_TUBES = wosc10.RES_WATER, wosc10.MATRIX_TUBES
SLOT_MALDI_TUBERACK = wosc10.SLOT_MALDI_TUBERACK
HS_SLOT, SLOT_SWAP = wosc10.HS_SLOT, wosc10.SLOT_SWAP
SLOT_TIPRACK_20, SLOT_DILUTION = wosc10.SLOT_TIPRACK_20, wosc10.SLOT_DILUTION
SLOT_TIPRACK_300, SLOT_RESERVOIR = wosc10.SLOT_TIPRACK_300, wosc10.SLOT_RESERVOIR
spread, MIX_RESERVOIR, MIX_LABEL = wosc10.spread, wosc10.MIX_RESERVOIR, wosc10.MIX_LABEL
ABTS_MIX_RESERVOIR = wosc10.ABTS_MIX_RESERVOIR
ABTS_MIX_LABEL = wosc10.ABTS_MIX_LABEL
RES_DEAD_UL = wosc10.RES_DEAD_UL

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import mm
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table,
                                 TableStyle, PageBreak, KeepTogether, HRFlowable,
                                 ListFlowable, ListItem)
from reportlab.lib.enums import TA_LEFT, TA_CENTER

# ---------------------------------------------------------------------------
# Build a concrete 3-enzyme worked example using the protocol's real functions,
# so every number in this guide is generated the same way the robot generates
# its own run log - nothing here is hand-typed or able to drift out of sync.
# ---------------------------------------------------------------------------
demo = [{'name': 'Lac-01', 'mg_ml': 2.40, 'mw_kda': 65.0},
        {'name': 'Lac-02', 'mg_ml': 1.85, 'mw_kda': 65.0},
        {'name': 'Lac-03', 'mg_ml': 3.10, 'mw_kda': 70.0}]
layout = build_layout(demo, to_um(2.40, 65.0))

# ---------------------------------------------------------------------------
# Styles
# ---------------------------------------------------------------------------
styles = getSampleStyleSheet()
NAVY = colors.HexColor('#1a2b4a')
GREY = colors.HexColor('#5a5a5a')
LIGHT = colors.HexColor('#eef1f6')
WARN_BG = colors.HexColor('#fdecea')
WARN_BORDER = colors.HexColor('#c62828')

styles.add(ParagraphStyle('DocTitle', parent=styles['Title'], fontSize=22,
                          textColor=NAVY, spaceAfter=4))
styles.add(ParagraphStyle('DocSubtitle', parent=styles['Normal'], fontSize=11,
                          textColor=GREY, spaceAfter=14))
styles.add(ParagraphStyle('H1', parent=styles['Heading1'], fontSize=14,
                          textColor=NAVY, spaceBefore=16, spaceAfter=6))
styles.add(ParagraphStyle('H2', parent=styles['Heading2'], fontSize=11.5,
                          textColor=NAVY, spaceBefore=10, spaceAfter=4))
styles.add(ParagraphStyle('Body', parent=styles['Normal'], fontSize=9.3,
                          leading=13, spaceAfter=6))
styles.add(ParagraphStyle('Small', parent=styles['Normal'], fontSize=8, leading=11,
                          textColor=GREY))
styles.add(ParagraphStyle('Cell', parent=styles['Normal'], fontSize=8.3, leading=10.5))
styles.add(ParagraphStyle('CellB', parent=styles['Cell'], fontName='Helvetica-Bold'))
styles.add(ParagraphStyle('WarnTitle', parent=styles['Normal'], fontSize=9.5,
                          fontName='Helvetica-Bold', textColor=WARN_BORDER))
styles.add(ParagraphStyle('MapMono', parent=styles['Normal'], fontName='Courier',
                          fontSize=7.6, leading=9.4))

PAGE_W, PAGE_H = A4
MARGIN = 12 * mm
CONTENT_W = PAGE_W - 2 * MARGIN


def p(text, style='Cell'):
    return Paragraph(text, styles[style])


def section(title):
    return [Paragraph(title, styles['H1']), HRFlowable(width='100%', thickness=1,
            color=NAVY, spaceAfter=8)]


def warn_box(title, lines):
    body = f'<b>{title}</b><br/>' + '<br/>'.join(lines)
    t = Table([[Paragraph(body, styles['Body'])]], colWidths=[CONTENT_W])
    t.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, -1), WARN_BG),
        ('BOX', (0, 0), (-1, -1), 0.75, WARN_BORDER),
        ('LEFTPADDING', (0, 0), (-1, -1), 10), ('RIGHTPADDING', (0, 0), (-1, -1), 10),
        ('TOPPADDING', (0, 0), (-1, -1), 8), ('BOTTOMPADDING', (0, 0), (-1, -1), 8),
    ]))
    return t


def styled_table(header, rows, col_widths, header_bg=NAVY, font_size=8.3,
                  row_bg=None):
    data = [[Paragraph(f'<b>{h}</b>', ParagraphStyle('h', parent=styles['Cell'],
             textColor=colors.white, fontSize=font_size)) for h in header]]
    for r in rows:
        data.append([c if hasattr(c, 'wrap') else
                    Paragraph(str(c), ParagraphStyle('c', parent=styles['Cell'],
                              fontSize=font_size)) for c in r])
    t = Table(data, colWidths=col_widths, repeatRows=1)
    style = [
        ('BACKGROUND', (0, 0), (-1, 0), header_bg),
        ('GRID', (0, 0), (-1, -1), 0.4, colors.HexColor('#c9c9c9')),
        ('VALIGN', (0, 0), (-1, -1), 'TOP'),
        ('LEFTPADDING', (0, 0), (-1, -1), 5), ('RIGHTPADDING', (0, 0), (-1, -1), 5),
        ('TOPPADDING', (0, 0), (-1, -1), 4), ('BOTTOMPADDING', (0, 0), (-1, -1), 4),
    ]
    if row_bg:
        for i in range(1, len(data)):
            if row_bg(i - 1):
                style.append(('BACKGROUND', (0, i), (-1, i), LIGHT))
    t.setStyle(TableStyle(style))
    return t


def swatch(hexcol):
    return Table([['']], colWidths=[9], rowHeights=[9],
                style=TableStyle([('BACKGROUND', (0, 0), (-1, -1),
                                   colors.HexColor(hexcol)),
                                  ('BOX', (0, 0), (-1, -1), 0.5, colors.black)]))


story = []

# ===========================================================================
# TITLE
# ===========================================================================
story.append(Paragraph('W-O-SC-10 Operator Setup Guide', styles['DocTitle']))
story.append(Paragraph(
    'NNBT/Purpald + guaiacol + ABTS + MALDI &mdash; laccase activity assay on the OT-2',
    styles['DocSubtitle']))
story.append(Paragraph(
    'This guide is a deck-setup reference for the person loading the robot. It does '
    'not replace the run log the protocol prints in the Opentrons app at the start of '
    'every run &ndash; that log has the <b>exact</b> volumes, wells and tube contents '
    'for the parameters chosen that day. This document explains the fixed structure '
    'that never changes: which slot holds what, which tube each reagent goes in, and '
    'the order operator actions happen in.', styles['Body']))

story.append(warn_box('Before you touch the deck', [
    'Open the protocol in the Opentrons app and fill in the runtime parameters '
    '(number of enzymes, Bradford concentrations, molar masses, tip start positions, '
    'MALDI on/off). <b>Do not edit the .py file.</b>',
    'The app\'s Setup screen will show a coloured, named liquid in every tube, '
    'reservoir well and plate well once parameters are set &ndash; use it to '
    'double-check this guide against your specific run before starting.',
]))

# ===========================================================================
# 1. OVERVIEW
# ===========================================================================
story += section('1. What this run does')
story.append(Paragraph(
    'One dilution series feeds three readouts from the same enzymes:', styles['Body']))
story.append(ListFlowable([
    ListItem(Paragraph('<b>NNBT plate</b> (slot 1, on the Heater-Shaker) &ndash; laccase '
             'activity via NNBT/Purpald, run in two parallel arms: plain and with '
             'guaiacol.', styles['Body'])),
    ListItem(Paragraph('<b>ABTS plate</b> (swapped into slot 5) &ndash; a fast kinetic '
             'activity readout at A414/A734, read immediately after spiking.',
             styles['Body'])),
    ListItem(Paragraph('<b>MALDI target</b> (swapped into slot 5 later) &ndash; optional. '
             'Spotted at t = '
             + ', '.join(f'{t:g}' for t in wosc10.MALDI_TIMEPOINTS_MIN)
             + ' min while the NNBT plate incubates. <b>One row per sample</b> '
             '(NC2, then the enzymes) with time running across the columns, two '
             'columns per timepoint (plain, then guaiacol).', styles['Body'])),
    ListItem(Paragraph('<b>10&times; readout plate</b> (swapped into slot 5 last) &ndash; '
             'an empty plate of the <b>same type as the NNBT plate</b>. Every finished '
             'reaction is diluted 10&times; into it after the incubation, and it '
             'finishes on the Heater-Shaker, whose adapter will not take a Corning '
             'flat plate.', styles['Body'])),
], bulletType='bullet', leftIndent=14))
story.append(Paragraph(
    f'The NNBT plate holds up to {wosc10.MAX_ENZYMES} enzymes plus two negative '
    'controls (NC2 = heat-inactivated enzyme, NC3 = active enzyme with no NNBT in the '
    f'mix) and a {len(wosc10.LAC_GRADIENT_UM)}-point lactaldehyde standard curve, each '
    'in triplicate, in both the plain and guaiacol arms. <b>NC1 (buffer only) is no '
    'longer on this plate</b> &ndash; ten enzymes need the row, and the 2026-09-15 run '
    'showed NC1 and NC2 agreeing to within one plate-noise SD, so NC2 (which is NC1 '
    'plus protein) is the better blank and does NC1\'s job. NC1 still runs on the '
    'ABTS plate, out of dilution well '
    f'{wosc10.ROW_LETTERS[len(PADA_NAMES)]}3.', styles['Body']))
story.append(Paragraph(
    '<b>The NNBT plate is read twice.</b> After the incubation, 15 &micro;L of every '
    'well goes into 135 &micro;L of milliQ on the readout plate. Both plates then get '
    'Purpald and develop on the Heater-Shaker, one after the other, and both are read '
    'at A530: save them as <b>"-undiluted"</b> and <b>"-10x"</b>. The 10&times; plate '
    'is the one whose standards should be on scale &ndash; on 2026-09-15 every '
    'standard from 500 to 6000 &micro;M read as one flat line on the blank, because '
    'the Purpald response peaks below 500 &micro;M and turns brown above it.',
    styles['Body']))

# ===========================================================================
# 2. DECK LAYOUT
# ===========================================================================
story += section('2. Deck layout')
deck_rows = [
    ['1', 'Heater-Shaker + NNBT plate (on the flat-bottom adapter)',
     'Stays on the deck the whole run. North (slot 4) must stay empty.'],
    ['2, 4, 9, 10, 11', 'EMPTY', 'Clearance, not spare space &ndash; leave empty.'],
    [SLOT_SWAP, 'Tube rack &rarr; ABTS plate &rarr; MALDI target (in that order)',
     'Only one of the three is ever on the deck at a time. The robot pauses and '
     'asks you to make each swap.'],
    [SLOT_TIPRACK_20, '20 &micro;L tip rack', '<b>p20 ONLY.</b> Never place a tube rack, '
     'reservoir or anything the p300 touches here (see &sect;7).'],
    [SLOT_DILUTION, 'Dilution plate (96-well)', 'Holds every diluted stock in columns '
     '1&ndash;3. Stays on the deck the whole run: the MALDI 1:5 dilutions go in '
     'columns 4&ndash;12 of the same plate.'],
    [SLOT_TIPRACK_300, '300 &micro;L tip rack &rarr; tube rack (MALDI only) &rarr; '
     '300 &micro;L tip rack', 'If MALDI is on, the robot asks you to take the tip rack '
     'out and put the tube rack back here with the capped matrix tube, then to swap '
     'them back before Purpald.'],
    [SLOT_RESERVOIR, '12-well reagent reservoir', 'North (slot 11) must stay empty.'],
    ['12', 'Fixed trash', 'Do not move or cover.'],
]
story.append(styled_table(['Slot', 'Contents', 'Notes'], deck_rows,
             [16 * mm, 68 * mm, 90 * mm]))
story.append(Spacer(1, 6))
story.append(Paragraph(
    'Load order for slot 5: put the <b>tube rack</b> on before starting the run. '
    'Keep the ABTS plate and MALDI target staged off-deck (labelled) until the robot '
    'pauses and asks for them by name.', styles['Small']))

# ===========================================================================
# 3. TIP RACKS
# ===========================================================================
story += section('3. Tip racks')
story.append(Paragraph(
    'Two full racks: a 20 &micro;L rack in slot 3 (p20, right mount) and a 300 &micro;L '
    'rack in slot 7 (p300, left mount). In the app, set <i>"20 uL tip row/column"</i> '
    'and <i>"300 uL tip column"</i> to the first unused tip if you are reusing partial '
    'racks &ndash; otherwise leave at A1. The app tells you at run start whether the '
    'loaded racks have enough tips; if not, keep a spare full rack of each size on '
    'hand, because the run will pause mid-way to ask for one.', styles['Body']))

# ===========================================================================
# 4. TUBE RACK
# ===========================================================================
story += section('4. Tube rack (slot 5, step 1 of 3)')
story.append(Paragraph(
    'A 24-position 1.5 mL tube rack (4 rows A&ndash;D &times; 6 columns). Enzyme stocks '
    'fill column-first: enzyme 1 &rarr; A1, enzyme 2 &rarr; B1, enzyme 3 &rarr; C1, '
    'enzyme 4 &rarr; D1, enzyme 5 &rarr; A2, and so on. Two positions are fixed and never '
    'move, regardless of how many enzymes are used:', styles['Body']))

tube_rows = [
    [f'{TUBE_HEAT_INACT} (fixed)', 'NC2 &ndash; heat-inactivated enzyme stock', swatch(COLOUR['nc2'])],
    [f'{TUBE_LACTALDEHYDE} (fixed)', 'Lactaldehyde stock, 100,000 &micro;M (1 M diluted 1:10)', swatch(COLOUR['lac'])],
]
for name, tube in zip(PADA_NAMES, PADA_TUBES):
    tube_rows.append([f'{tube} (fixed)',
                      f'{name} &ndash; ABTS positive control, made up by hand from '
                      'powder',
                      swatch(COLOUR['pada'])])
tube_rows.append([', '.join(MATRIX_TUBES),
                  'MALDI matrix (DHB/acetonitrile), <b>CAPPED</b> &ndash; only if MALDI '
                  f'is on. Added when the rack goes back into slot {SLOT_MALDI_TUBERACK}; '
                  'the second tube only if the run log asks for it.',
                  swatch(COLOUR['matrix'])])
story.append(styled_table(['Tube', 'Contents', 'Colour'], tube_rows,
             [30 * mm, 130 * mm, 14 * mm]))
story.append(Spacer(1, 6))
story.append(Paragraph(
    'Enzyme positions (A1, B1, C1, D1, A2 &hellip;) hold each enzyme\'s own stock at '
    'whatever concentration your Bradford assay measured &ndash; the robot dilutes it '
    'to the common target concentration itself. The enzyme chosen as <i>"NC3 uses '
    'enzyme #"</i> is drawn from its own tube a second time, so nothing is short-pipetted.'
    ' <b>Exact tube volumes needed are in the run log</b> &ndash; they depend on how '
    'many enzymes you run and the tip use they need to survive.', styles['Body']))

# ===========================================================================
# 5. RESERVOIR
# ===========================================================================
story += section('5. Reagent reservoir (slot 8)')
story.append(Paragraph(
    'A 12-well reservoir. All reagent mixes below are <b>premixed by hand off-deck</b> '
    'before the run &ndash; the robot only dispenses them. Pour to the volume shown in '
    'the run log (it already includes dead volume); the table below shows which well '
    'gets which reagent.', styles['Body']))

res_rows = [
    (RES_NNBT, 'NNBT reaction mix (assay buffer + NNBT)', 'nnbt'),
    (RES_BUFFER, 'Assay / dilution buffer', 'buffer'),
    (RES_PURPALD, 'Purpald reagent', 'purpald'),
    (RES_ABTS, 'Laccase ABTS reaction mix', 'abts_mix'),
    (RES_PADA, 'PaDa-1 ABTS reaction mix &ndash; H<sub>2</sub>O<sub>2</sub> instead of '
     'Cu<sup>2+</sup>, and a different pH. <b>Not interchangeable with the laccase '
     'mix.</b>', 'pada_mix'),
    (RES_GUAIACOL, 'Guaiacol reaction mix (buffer + NNBT + guaiacol)', 'gua'),
    (RES_NO_NNBT, 'No-NNBT mix (NC3, plain arm)', 'no_nnbt'),
    (RES_NO_NNBT_GUA, 'No-NNBT + guaiacol mix (NC3, guaiacol arm)', 'no_nnbt'),
    (RES_WATER, 'milliQ water &ndash; the 10&times; readout diluent (always), plus the '
     'MALDI 1:5 premix if spotting is on. <b>Both wells</b>: one cannot hold it.',
     'water'),
]
res_table_rows = []
for wells, label, tag in res_rows:
    wells_list = wells if isinstance(wells, list) else [wells]
    res_table_rows.append([', '.join(wells_list), label, swatch(COLOUR[tag])])
story.append(styled_table(['Well(s)', 'Reagent', 'Colour'], res_table_rows,
             [24 * mm, 130 * mm, 20 * mm]))
story.append(Spacer(1, 4))
story.append(Paragraph(
    'Wells shown as two positions (e.g. NNBT mix, guaiacol mix) mean the reagent is '
    'split across both if a single 15 mL well would not hold enough for a large run '
    '&ndash; the run log tells you whether the second well is actually needed.',
    styles['Small']))

# ===========================================================================
# 6. RUN SEQUENCE
# ===========================================================================
story += section('6. Run sequence &amp; where the robot will pause')
story.append(Paragraph(
    'Steps 1&ndash;5, 7&ndash;9, 11, 13&ndash;14 and 16&ndash;17 run automatically. '
    'The robot stops and waits for you at every step marked <b>PAUSE</b> &ndash; read '
    'the on-screen message before resuming.', styles['Body']))

seq_rows = [
    ['1', 'Buffer dispensed into every dilution well', ''],
    ['2', 'Enzyme, NC3 and heat-inactivated stocks diluted and mixed to target conc.', ''],
    ['3', 'Lactaldehyde standard curve built by serial dilution', ''],
    ['4', 'PaDa-1 1:1000 moved from its tube into a dilution well',
     'rack is about to leave'],
    ['5', '10 &micro;L spikes into all four NNBT plate blocks', ''],
    ['6', 'PAUSE', '<b>Swap the tube rack for the ABTS plate in slot 5.</b>'],
    ['7', '10 &micro;L spikes into the ABTS plate, one sample at a time', ''],
    ['8', 'Reaction mixes dispensed into the NNBT plate', ''],
    ['9', 'Both ABTS mixes dispensed by column &ndash; always last, since ABTS reacts '
     'on contact', 'column 1 takes the PaDa-1 mix, the rest the laccase mix'],
    ['10', 'PAUSE', '<b>Take the ABTS plate to the reader immediately</b> (read '
     'A414/A734), then seal the NNBT plate. If MALDI is on: MALDI target into slot 5, '
     '300 &micro;L tips out of slot 7, tube rack with the capped matrix tube into slot 7.'],
    ['11', 'Incubate 2 h @ 40&deg;C', 'If MALDI is on, pauses at t = '
     + ', '.join(f'{t:g}' for t in wosc10.MALDI_TIMEPOINTS_MIN) + ' min: '
     '<b>unseal &rarr; open the matrix tube once for the whole session, close it at '
     'the end &rarr; reseal.</b> Not once per spot.'],
    ['12', 'PAUSE for the 10&times; readout',
     '<b>Remove the plate seal. Put the empty readout plate in slot 5 and a FRESH '
     '20 &micro;L tip rack in slot 3</b> &ndash; the 10&times; pass spends a whole '
     'rack, one tip column per plate column. Check the run log for how many racks '
     'this run needs in total.'],
    ['13', '135 &micro;L milliQ into the readout plate, then 15 &micro;L out of every '
     'NNBT well into it, mixed', ''],
    ['14', '45 &micro;L Purpald into the <b>undiluted</b> plate, developed 10 min',
     'this plate is quenched first on purpose &ndash; it still holds active enzyme at '
     'full NNBT, while the 10&times; plate already runs ten times slower'],
    ['15', 'PAUSE to swap the plates',
     '<b>Take the undiluted plate off the Heater-Shaker to the reader (A530, save as '
     '"-undiluted"), then put the 10&times; plate onto the Heater-Shaker.</b>'],
    ['16', '10&times; plate shaken 2 min at 1000 rpm', '15 &micro;L under 135 &micro;L '
     'does not mix itself, and Purpald meeting a concentrated bolus browns it for good'],
    ['17', '50 &micro;L Purpald into the 10&times; plate, developed 10 min, done',
     '<b>Read A530, save as "-10x".</b>'],
]
story.append(styled_table(['#', 'Action', 'Operator note'], seq_rows,
             [10 * mm, 78 * mm, 96 * mm],
             row_bg=lambda i: seq_rows[i][0] in ('6', '10', '12', '15')))

story.append(PageBreak())

# ===========================================================================
# 7. SAFETY RULES
# ===========================================================================
story += section('7. Three rules that have already cost a run')
story.append(warn_box('Rule 1 &ndash; Slot 9 can never take a single-nozzle move', [
    'The fixed trash in slot 12 sits directly north of slot 9, and the simulator does '
    'not model that clearance. Never place labware needing single-tip access in slot 9. '
    'The protocol checks this automatically and will refuse to load if violated.']))
story.append(Spacer(1, 6))
story.append(warn_box('Rule 2 &ndash; a HEIGHT is not a DEPTH', [
    'A tip positioned above the liquid surface draws air; a tip pressed to a flat well '
    'bottom seals against it and moves nothing. Every liquid-handling height in this '
    'protocol is computed from the well\'s remaining volume &ndash; if you ever '
    'manually jog the robot, do not assume a fixed height works at every fill level.']))
story.append(Spacer(1, 6))
story.append(warn_box('Rule 3 &ndash; the LEFT mount (p300) cannot reach the right-hand deck column', [
    'Both pipettes ride the same X carriage ~34 mm apart, against a hard limit at '
    'X 418 mm, so the left-mount p300 tops out near deck X 382. Slots 3/6/9 start at '
    'x = 265, so anything past x &asymp; 117 within one of those slots is physically '
    'unreachable by the p300 &ndash; and a touch-tip adds another 4 mm on top.',
    '<b>This is what caused the 2026-09-11 run failure:</b> the tube rack was in '
    'slot 3 and a touch-tip on tube D6 commanded X422.2, tripping a hard-limit alarm.',
    '<b>Operator rule: slot 3 must only ever hold the 20 &micro;L tip rack.</b> Never '
    'place a tube rack, reservoir, or any labware the p300 touches there. The tube '
    'rack, ABTS plate and MALDI target all belong in slot 5, which sits in the middle '
    'deck column and is reachable by both pipettes. During MALDI only, the tube rack '
    'goes in slot 7, where only the p20 visits it.']))

# ===========================================================================
# 8. WORKED EXAMPLE
# ===========================================================================
story += section('8. Worked example (3 enzymes) &ndash; for illustration only')
story.append(Paragraph(
    'The tables below are generated from the protocol\'s own layout code using 3 '
    'placeholder enzymes, purely to show the <i>shape</i> of the plate maps and what a '
    'completed reservoir/tube-rack setup looks like. <b>Do not pour these volumes for '
    'a real run</b> &ndash; use the run log the app prints for your actual parameters.',
    styles['Body']))

story.append(Paragraph('NNBT plate map (slot 1)', styles['H2']))
story.append(Paragraph(
    'N = NNBT mix &middot; G = guaiacol mix &middot; x = no-NNBT mix (NC3) &middot; '
    'y = no-NNBT + guaiacol mix (NC3). Number = group index, keyed to the table below.',
    styles['Small']))

pos = {w: i + 1 for i, g in enumerate(layout['nnbt']) for w in g['wells']}
mix = {w: g['mix'] for g in layout['nnbt'] for w in g['wells']}
SYMBOL = {'nnbt': 'N', 'gua': 'G', 'no_nnbt': 'x', 'no_nnbt_gua': 'y'}
ROWS = list('ABCDEFGH')
header = [''] + [str(c) for c in range(1, 13)]
map_rows = []
for r in ROWS:
    row = [r]
    for c in range(1, 13):
        w = f'{r}{c}'
        row.append(f'{SYMBOL[mix[w]]}{pos[w]:02d}' if w in pos else '.')
    map_rows.append(row)
story.append(styled_table(header, map_rows,
             [9 * mm] + [13.6 * mm] * 12, font_size=7.6))
story.append(Spacer(1, 6))

nnbt_group_rows = [[str(i + 1), ','.join(g['wells']), g['label']]
                   for i, g in enumerate(layout['nnbt'])]
story.append(styled_table(['#', 'Wells', 'Content'], nnbt_group_rows,
             [10 * mm, 40 * mm, 134 * mm], font_size=7.8))

story.append(PageBreak())
story.append(Paragraph('ABTS plate map (slot 5, swapped in after the dilutions)',
             styles['H2']))
abts_rows = [[str(i + 1), ','.join(g['wells']), g['label'],
              ABTS_MIX_LABEL[g['mix']]]
            for i, g in enumerate(layout['abts'])]
story.append(styled_table(['#', 'Wells', 'Content', 'Reaction mix'], abts_rows,
             [10 * mm, 34 * mm, 90 * mm, 50 * mm], font_size=7.8))
story.append(Spacer(1, 6))
story.append(Paragraph(
    'The two PaDa-1 strengths are the positive control &ndash; there is no ABTS radical '
    'standard curve. They run <i>down</i> column 1 because the reaction mixes go in by '
    'whole columns, and a column can only take one mix. Wells in a part-filled column '
    'still receive reaction mix but hold no sample: <b>they are blank, do not read '
    'them.</b>', styles['Body']))

story.append(Spacer(1, 10))
story.append(Paragraph('Dilution plate (slot 6), 150 &micro;L per well', styles['H2']))
dil_rows = []
for d in layout['dilutions']:
    origin = {'tube': 'tube ', 'res': 'res '}.get(d['src_kind'], 'dil ') + d['src']
    dil_rows.append([d['well'], origin, f'{d["stock_vol"]:.1f}', f'{d["buffer_vol"]:.1f}',
                     d['name']])
story.append(styled_table(['Well', 'From', 'Stock µL', 'Buffer µL', 'Content'],
             dil_rows, [12 * mm, 24 * mm, 18 * mm, 18 * mm, 112 * mm], font_size=7.6))

story.append(Spacer(1, 10))
story.append(Paragraph('Reservoir (slot 8) &ndash; example pour volumes', styles['H2']))
# Same computation render() uses internally - not text-scraped, so it can't pick up
# a false match from another table (e.g. the dilution plate also has a well "A1").
needs = {RES_BUFFER: ('assay buffer', layout['buffer_total'])}
for key, total in layout['per_mix'].items():
    for w, v in spread(total, MIX_RESERVOIR[key], MIX_LABEL[key]).items():
        needs[w] = (MIX_LABEL[key], v)
for w, v in spread(layout['purpald_total'], wosc10.RES_PURPALD, 'Purpald').items():
    needs[w] = ('Purpald reagent', v)
for key, total in layout['abts_per_mix'].items():
    for w, v in spread(total, ABTS_MIX_RESERVOIR[key], ABTS_MIX_LABEL[key]).items():
        needs[w] = (ABTS_MIX_LABEL[key], v)
res_example_rows = []
for w in sorted(needs, key=lambda x: int(x[1:])):
    name, vol = needs[w]
    res_example_rows.append([w, name, f'{(vol + RES_DEAD_UL) / 1000:.1f} mL'])
story.append(styled_table(['Well', 'Reagent', 'Volume'], res_example_rows,
             [16 * mm, 100 * mm, 68 * mm], font_size=8))

story.append(Spacer(1, 10))
story.append(Paragraph('Tube rack &ndash; example volumes needed', styles['H2']))
tubes = {}
for d in layout['dilutions']:
    if d['src_kind'] != 'tube':
        continue
    vol, names = tubes.get(d['src'], (0.0, []))
    tubes[d['src']] = (vol + (d['stock_vol'] or 40.0), names + [d['name']])
tube_example_rows = []
for tube in sorted(tubes, key=lambda t: (int(t[1:]), t[0])):
    vol, names = tubes[tube]
    tube_example_rows.append([tube, ' + '.join(names), f'≥ {vol + 300:.0f} µL'])
story.append(styled_table(['Tube', 'Contents', 'Min. volume'], tube_example_rows,
             [16 * mm, 118 * mm, 50 * mm], font_size=8))

story.append(Spacer(1, 14))
story.append(Paragraph(
    'Generated from Lac_W-O-SC-10_nnbt_guaiacol_abts_maldi_ot2.py &ndash; '
    'Kristine Toft Johansen (s215098). This guide reflects the protocol as of the '
    '"slot-swap left-mount reach" fix (tube rack confined to slot 5; slot 3 reserved '
    'for the 20 &micro;L tip rack only).', styles['Small']))

# ===========================================================================
doc = SimpleDocTemplate('WOSC10_Operator_Setup_Guide.pdf', pagesize=A4,
                        topMargin=14 * mm, bottomMargin=14 * mm,
                        leftMargin=MARGIN, rightMargin=MARGIN,
                        title='W-O-SC-10 Operator Setup Guide')
doc.build(story)
print('Wrote WOSC10_Operator_Setup_Guide.pdf')
