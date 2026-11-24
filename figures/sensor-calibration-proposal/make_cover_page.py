"""Cover / rationale pages for the proposal PDF."""
import os
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'DejaVu Sans']
MM = 1 / 25.4
OUT = os.path.dirname(os.path.abspath(__file__))

TITLE = 'Proposed new analyses & visualizations'
SUB = ('Paragraph: "A sensor-informed analytical framework for robust calibration and\n'
       'identification of context-specific drivers and dependencies"')

BODY = [
    ('THE GAP', 'b', [
        'The paragraph makes four quantitative claims that the current panels assert but never show:',
        '  1. Sensors "allow quantitative calibration" — but no panel shows calibration changing an answer.',
        '  2. The 20% cutoff is "a dynamic functional threshold akin to an analog lever" — shown at one',
        '      value (Fig. S9D), so the reader cannot see the lever actually move.',
        '  3. Editing is "highly correlated across conditions" — shown as distributions, not per-gRNA agreement.',
        '  4. Larger library ⇒ "decrease in biological signal relative to noise" — stated, never measured.',
        '',
        'Every panel below is computed from your existing tables (LFC-FDR/*_OG_LFC_FDR.csv, counts/og/*).',
        'No new experiments. Numbers reproduce your text: 803/1,034 ABE gRNAs pass ≥20% target editing at D15.',
    ]),
    ('MAIN FIGURE — calibration as a tunable lever', 'b', [
        'a  Editing-outcome composition (perfect / bystander / unedited) per editor × tissue. Replaces the',
        '    two-bar Fig. 2D with a composition that shows what the other ~40% of reads actually are.',
        'b  Editing purity (perfect ÷ target) — ABE median ~60%, CBE ~11%. A single number justifying',
        '    why you report two metrics at all.',
        'c  Per-gRNA editing reproducibility ACROSS tissues (ρ = 0.85–0.96). This is the real evidence for',
        '    "platform robustness" — same gRNA, same efficiency, different mouse and organ.',
        'd  ** THE KEY NEW PANEL ** Signal-to-noise (SD targeting ÷ SD control gRNAs) vs editing cutoff.',
        '    Rises monotonically 0→90% in every ABE context. This IS the analog lever, drawn.',
        'e  Inter-replicate ρ vs cutoff. Rises for ABE (BM 0.41→0.51, D15 0.74→0.86); FLAT for CBE.',
        'f  The cost of the lever: gRNAs retained vs cutoff — the sensitivity/stringency trade-off.',
        'g  LFC vs editing with the excluded zone shaded: shows calibration changing the hit list.',
        'h  Effective library coverage (Shannon) — bottlenecking as a % of library actually surviving.',
        'i  Coverage vs reproducibility: one scatter tying library size to recoverable signal.',
    ]),
    ('EXTENDED FIGURES', 'b', [
        'E1  Threshold-sweep atlas (ABE top, CBE bottom) × SNR / % significant / % retained.',
        '     Headline: ABE SNR > 1 at every cutoff and climbing; CBE SNR < 1 at EVERY cutoff —',
        '     no sensor filter can rescue signal that bottlenecking already destroyed. This is a much',
        '     sharper claim than "a concomitant decrease in biological signal relative to noise".',
        'E2  Target vs perfect editing: hexbin agreement, purity distributions, and SNR under each',
        '     metric. ABE purity is strongly BIMODAL — gRNAs are near-pure or heavily contaminated,',
        '     few in between. Argues for treating purity as a class, not a continuum.',
        'E3  Bottlenecking: Lorenz curves, Gini, effective gRNA number (log), skew vs reproducibility,',
        '     control-gRNA noise floor, and dropout. Panel d spans two orders of magnitude of skew.',
        'E4  Specificity control: editing depth and sensor read depth do NOT manufacture effect size',
        '     (ρ ≈ −0.2 to −0.3, i.e. if anything the opposite). Pre-empts the obvious reviewer',
        '     objection that the cutoff is self-fulfilling.',
    ]),
    ('CAVEATS — read before using', 'b', [
        '• DRAFT panels for discussion, not final figures. Styling follows your sci-plot conventions',
        '   (7pt, open axes, fixed tissue/editor colors) but panel order and composition are proposals.',
        '• SNR is defined as SD(targeting LFC) ÷ SD(control LFC). Controls = safe- + non-targeting.',
        '   ABE has only 36 control gRNAs, so the ABE null is estimated from few guides — worth a',
        '   bootstrap CI before this goes in a paper.',
        '• CBE SNR < 1 is a real feature of the data, but check it is not driven by the control gRNAs',
        '   themselves being bottlenecked (they are in the same pool). E3e shows the control noise',
        '   floor is genuinely wider in CBE, which supports the interpretation.',
        '• Panels use sensor reads ≥ 100 and FDR (Stouffer). Meninges ABE n is small in places.',
        '• Fig. S9A (position-wise +5..+7 editing window) is NOT reproduced here — that needs the',
        '   per-position allele tables, not the compact CRISPResso summaries. Happy to add it.',
    ]),
]


def cover_fig():
    blocks = BODY
    fig = plt.figure(figsize=(183 * MM, 247 * MM))
    fig.text(.02, .975, TITLE, fontsize=13, fontweight='bold', va='top')
    fig.text(.02, .952, SUB, fontsize=7, va='top', color='#333333')
    fig.text(.02, .916, 'B-ALL in vivo base editing · draft for review · not final figures',
             fontsize=6.5, va='top', color='#777777', style='italic')
    y = .885
    for head, _, lines in blocks:
        fig.text(.02, y, head, fontsize=8.5, fontweight='bold', va='top',
                 color='#1a3c6e')
        y -= .021
        for ln in lines:
            fig.text(.02, y, ln, fontsize=6.6, va='top', family='DejaVu Sans')
            y -= .0152
        y -= .012
    return fig


if __name__ == '__main__':
    f = cover_fig()
    f.savefig(f'{OUT}/cover.pdf', bbox_inches='tight')
    print('cover done')
