"""
The two figures Joy asked to have fixed on the 7 November call.

ASK 2, the dot figure. The 4 November version drew a maternal loss and a paternal loss as two
identical clouds and told the reader that red dots were the only thing separating them. On a raw
B-allele-frequency axis that is true, and that is the defect: the figure never draws the axis the
program actually uses. It also coloured only one of the two informative marker classes, so the
other class sat in the opposite pile drawn identically in black.

ASK 1, the share scale. The same slide named the quantity "B allele frequency, averaged over about
a thousand markers". It is the paternal share, it is a median, and it is judged against the
sample's own centre rather than against 0.50. Its two gain labels were also swapped: gainOrigin
asserts that a share of 0.667 is a PATERNAL gain and 0.333 a maternal one.

Both figures are schematics of the model, not plots of an array. Informative markers only, which is
why a normal embryo sits in one pile rather than three.

Run: python3 tools/joy_figures.py [outdir]
"""
import sys, os, random
from reportlab.pdfgen import canvas
from reportlab.lib.colors import Color

W, H = 792, 612
INK   = Color(0.10, 0.10, 0.11)
MUTE  = Color(0.42, 0.42, 0.45)
RULE  = Color(0.80, 0.80, 0.82)
# The two informative marker classes. Distinguished by shape as well as colour so the figure
# survives being printed in grey, which is the project's own rule for its product figures.
CLS_A = Color(0.13, 0.31, 0.56)   # father AA, mother BB   filled
CLS_B = Color(0.72, 0.34, 0.10)   # father BB, mother AA   hollow

random.seed(20261108)


def text(c, x, y, s, size=9, colour=INK, font="Helvetica", anchor="start"):
    c.setFillColor(colour); c.setFont(font, size)
    if anchor == "middle": c.drawCentredString(x, y, s)
    elif anchor == "end":  c.drawRightString(x, y, s)
    else:                  c.drawString(x, y, s)


def panel(c, x0, y0, w, h, title, points, xlabel, xticks, ylo=-1.35, yhi=0.45):
    """One scatter panel. `points` is a list of (xval, logr, cls)."""
    c.setStrokeColor(RULE); c.setLineWidth(0.7)
    c.rect(x0, y0, w, h, stroke=1, fill=0)
    text(c, x0 + w / 2, y0 + h + 9, title, 9.5, INK, "Helvetica-Bold", "middle")

    sx = lambda v: x0 + 14 + v * (w - 28)
    sy = lambda v: y0 + 10 + (v - ylo) / (yhi - ylo) * (h - 20)

    # zero line, the sample's own median
    c.setStrokeColor(RULE); c.setDash(2, 2); c.setLineWidth(0.6)
    c.line(x0 + 4, sy(0), x0 + w - 4, sy(0)); c.setDash()

    for xv, lr, cls in points:
        col = CLS_A if cls == "A" else CLS_B
        c.setStrokeColor(col); c.setLineWidth(0.8)
        if cls == "A":
            c.setFillColor(col); c.circle(sx(xv), sy(lr), 2.5, stroke=0, fill=1)
        else:
            c.circle(sx(xv), sy(lr), 2.5, stroke=1, fill=0)

    c.setStrokeColor(MUTE); c.setLineWidth(0.6)
    for t, lab in xticks:
        c.line(sx(t), y0, sx(t), y0 - 3)
        text(c, sx(t), y0 - 12, lab, 7.5, MUTE, "Helvetica", "middle")
    text(c, x0 + w / 2, y0 - 24, xlabel, 8, MUTE, "Helvetica", "middle")
    for t in (0.0, -1.0):
        text(c, x0 - 4, sy(t) - 2.5, f"{t:+.0f}" if t else "0", 7.5, MUTE, "Helvetica", "end")


def cloud(xa, xb, logr, n=34, jx=0.018, jy=0.055):
    """n markers, half of each class, at the given x for class A and class B."""
    pts = []
    for i in range(n):
        cls = "A" if i % 2 == 0 else "B"
        x = (xa if cls == "A" else xb) + random.gauss(0, jx)
        pts.append((min(max(x, 0.004), 0.996), logr + random.gauss(0, jy), cls))
    return pts


def fig_dots(path):
    c = canvas.Canvas(path, pagesize=(W, H))
    text(c, W / 2, H - 42, "Two different losses, and the axis that tells them apart",
         15, INK, "Helvetica-Bold", "middle")
    text(c, W / 2, H - 60,
         "Informative markers only: the father and the mother are homozygous for different "
         "alleles, so a normal embryo must read AB.", 9, MUTE, "Helvetica", "middle")

    cols = [("Both copies there", 0.5, 0.5, 0.0),
            ("The mother's copy is gone", 0.0, 1.0, -1.0),
            ("The father's copy is gone", 1.0, 0.0, -1.0)]
    pw, gap, left = 215, 24, 59
    xt_baf = [(0.0, "0"), (0.5, "0.5"), (1.0, "1")]

    # ROW 1, what the array reports
    text(c, left, 494, "WHAT THE ARRAY REPORTS", 9, INK, "Helvetica-Bold")
    text(c, left + 142, 494, "the two losses are the same picture", 9, MUTE)
    for i, (t, xa, xb, lr) in enumerate(cols):
        panel(c, left + i * (pw + gap), 350, pw, 108, t, cloud(xa, xb, lr),
              "B allele frequency", xt_baf)

    # ROW 2, what the program reads
    text(c, left, 284, "WHAT THE PROGRAM READS", 9, INK, "Helvetica-Bold")
    text(c, left + 146, 284, "one flip per marker, and they separate completely", 9, MUTE)
    for i, (t, xa, xb, lr) in enumerate(cols):
        # the flip: class A keeps 1 - BAF, class B keeps BAF. Both classes then put the
        # father's allele at 1.
        pts = [(1 - x if cls == "A" else x, y, cls) for x, y, cls in cloud(xa, xb, lr)]
        panel(c, left + i * (pw + gap), 152, pw, 108, "", pts,
              "paternal share", [(0.0, "0"), (0.5, "0.5"), (1.0, "1")])

    c.setStrokeColor(RULE); c.setLineWidth(0.7); c.line(left, 112, W - left, 112)

    # legend
    c.setFillColor(CLS_A); c.circle(left + 4, 92, 2.5, stroke=0, fill=1)
    text(c, left + 14, 89, "father AA, mother BB", 8.5, INK)
    c.setStrokeColor(CLS_B); c.setLineWidth(0.8); c.circle(left + 150, 92, 2.5, stroke=1, fill=0)
    text(c, left + 160, 89, "father BB, mother AA", 8.5, INK)

    text(c, left, 68,
         "The flip uses the parents' genotypes, one marker at a time: where the father is AA the "
         "program keeps 1 - BAF, where he is BB it keeps BAF.", 9, INK)
    text(c, left, 54,
         "Both arrangements then put the father's allele at 1. A lost maternal copy leaves his "
         "allele alone and the share goes to 1; a lost paternal copy sends it to 0.", 9, INK)
    text(c, left, 36,
         "No colour is needed to tell the two apart. The colours above only show that the same "
         "two piles are built from opposite marker classes.", 9, MUTE)
    c.showPage(); c.save()


def fig_scale(path):
    c = canvas.Canvas(path, pagesize=(W, H))
    text(c, W / 2, H - 44, "What the share is, and what its values mean",
         15, INK, "Helvetica-Bold", "middle")
    text(c, W / 2, H - 62,
         "It is not a B allele frequency. It is a B allele frequency flipped per marker by the "
         "parents' genotypes, then taken as a median.", 9, MUTE, "Helvetica", "middle")

    # the definition, as three steps
    x0, y = 70, 468
    for n, head, body in [
        ("1", "Keep only the markers that can speak",
         "The father and the mother are homozygous for different alleles. 11 markers per Mb, "
         "about 1,000 per chromosome."),
        ("2", "Flip each marker so high always means father",
         "Father AA, mother BB: keep 1 - BAF.      Father BB, mother AA: keep BAF."),
        ("3", "Take the median over the region, and compare it to this sample's own centre",
         "A median, not a mean: the distribution has a homozygous tail that a mean chases. The "
         "centre is the sample's own genome-wide median, not 0.50."),
    ]:
        c.setFillColor(INK); c.circle(x0 + 5, y + 3, 8, stroke=0, fill=1)
        text(c, x0 + 5, y, n, 9, Color(1, 1, 1), "Helvetica-Bold", "middle")
        text(c, x0 + 22, y, head, 10.5, INK, "Helvetica-Bold")
        text(c, x0 + 22, y - 14, body, 9, MUTE)
        y -= 46

    # the scale
    bx0, bx1, by = 110, W - 110, 250
    c.setStrokeColor(INK); c.setLineWidth(1.0); c.line(bx0, by, bx1, by)
    sx = lambda v: bx0 + v * (bx1 - bx0)
    marks = [(0.0, "0", "the father's copy is gone", "his alleles are not there at all"),
             (1/3, "0.33", "an extra MATERNAL copy", "two of hers against one of his"),
             (0.5, "0.50", "one from each parent", "every informative marker reads AB"),
             (2/3, "0.67", "an extra PATERNAL copy", "two of his against one of hers"),
             (1.0, "1.00", "the mother's copy is gone", "his allele stands alone")]
    for v, lab, head, sub in marks:
        up = v in (1/3, 2/3)
        c.setStrokeColor(INK); c.setLineWidth(1.0)
        c.line(sx(v), by - 5, sx(v), by + 5)
        c.setFillColor(INK); c.circle(sx(v), by, 3.2, stroke=0, fill=1)
        text(c, sx(v), by - 20, lab, 11, INK, "Helvetica-Bold", "middle")
        ty = by + 46 if up else by - 40
        text(c, sx(v), ty, head, 8.6, INK, "Helvetica-Bold", "middle")
        text(c, sx(v), ty - 11, sub, 8, MUTE, "Helvetica", "middle")
        c.setStrokeColor(RULE); c.setLineWidth(0.6); c.setDash(2, 2)
        c.line(sx(v), by + (10 if up else -26), sx(v), ty - (22 if up else -16))
        c.setDash()
    text(c, W / 2, by - 72, "paternal share", 9, MUTE, "Helvetica", "middle")

    c.setStrokeColor(RULE); c.setLineWidth(0.7); c.line(70, 132, W - 70, 132)
    text(c, 70, 110, "A loss is the bigger signal of the two.", 9.5, INK, "Helvetica-Bold")
    text(c, 70, 96,
         "Losing one of two copies removes that parent's alleles entirely and sends the share to "
         "0 or 1. Adding a third copy only moves it to 0.33 or 0.67,", 9, INK)
    text(c, 70, 82,
         "so a gain is read against a margin rather than against a pole. Amplification compresses "
         "the observed value, and a partial or mosaic event lands between the two.", 9, INK)
    text(c, 70, 58,
         "Source: paternalShare and callGainOrigin in web/src/gainOrigin.ts, and the directions "
         "asserted in gainOrigin.check.ts.", 8.5, MUTE)
    c.showPage(); c.save()


def wrap(c, x, y, w, s, size=8.6, colour=MUTE, lead=11, font="Helvetica"):
    """Greedy wrap. The boxes are narrow and the text is written to fit, not reflowed by hand."""
    c.setFont(font, size)
    line, out = "", []
    for word in s.split():
        t = (line + " " + word).strip()
        if c.stringWidth(t, font, size) > w and line:
            out.append(line); line = word
        else:
            line = t
    if line: out.append(line)
    for i, l in enumerate(out):
        text(c, x, y - i * lead, l, size, colour, font)
    return y - (len(out) - 1) * lead


def box(c, x, y, w, h, fill=None, edge=RULE, lw=0.8):
    if fill: c.setFillColor(fill); c.rect(x, y, w, h, stroke=0, fill=1)
    c.setStrokeColor(edge); c.setLineWidth(lw); c.rect(x, y, w, h, stroke=1, fill=0)


def arrow(c, x0, x1, y):
    c.setStrokeColor(MUTE); c.setLineWidth(1.1)
    c.line(x0, y, x1 - 5, y)
    c.setFillColor(MUTE)
    pth = c.beginPath(); pth.moveTo(x1, y); pth.lineTo(x1 - 7, y + 3.4)
    pth.lineTo(x1 - 7, y - 3.4); pth.close(); c.drawPath(pth, stroke=0, fill=1)


def fig_schematic(path):
    c = canvas.Canvas(path, pagesize=(W, H))
    text(c, W / 2, H - 48, "What the program does with one array",
         16, INK, "Helvetica-Bold", "middle")
    text(c, W / 2, H - 68, "Three steps. Every threshold lives in the code, not in this picture.",
         9.5, MUTE, "Helvetica", "middle")

    bw, by, bh = 214, 276, 190
    xs = [52, 52 + bw + 32, 52 + 2 * (bw + 32)]
    top = by + bh

    # 1
    box(c, xs[0], by, bw, bh)
    text(c, xs[0] + 14, top - 24, "1", 13, INK, "Helvetica-Bold")
    text(c, xs[0] + 30, top - 24, "Pick the markers", 11, INK, "Helvetica-Bold")
    text(c, xs[0] + 30, top - 38, "that can speak", 11, INK, "Helvetica-Bold")
    wrap(c, xs[0] + 14, top - 62, bw - 28,
         "The father and the mother are homozygous for different alleles, so a normal embryo "
         "has to read AB.", colour=INK)
    wrap(c, xs[0] + 14, by + 46, bw - 28,
         "Every other marker reads the same whichever parent a copy came from, so it says "
         "nothing about parent of origin.")

    # 2, two lanes with the divider placed between them rather than through them
    box(c, xs[1], by, bw, bh)
    text(c, xs[1] + 14, top - 24, "2", 13, INK, "Helvetica-Bold")
    text(c, xs[1] + 30, top - 24, "Read two things", 11, INK, "Helvetica-Bold")
    text(c, xs[1] + 30, top - 38, "at those markers", 11, INK, "Helvetica-Bold")
    text(c, xs[1] + 14, top - 62, "WHICH PARENT", 8.5, CLS_A, "Helvetica-Bold")
    wrap(c, xs[1] + 14, top - 76, bw - 28,
         "The share of the signal that is the father's.", colour=INK)
    c.setStrokeColor(RULE); c.setLineWidth(0.6)
    c.line(xs[1] + 14, top - 98, xs[1] + bw - 14, top - 98)
    text(c, xs[1] + 14, top - 116, "ADDED OR LOST", 8.5, CLS_B, "Helvetica-Bold")
    wrap(c, xs[1] + 14, top - 130, bw - 28,
         "The intensity, against this sample's own middle.", colour=INK)
    wrap(c, xs[1] + 14, by + 34, bw - 28,
         "Neither on its own is enough. Intensity cannot name a parent.")

    # 3
    box(c, xs[2], by, bw, bh)
    text(c, xs[2] + 14, top - 24, "3", 13, INK, "Helvetica-Bold")
    text(c, xs[2] + 30, top - 24, "Put them together", 11, INK, "Helvetica-Bold")
    wrap(c, xs[2] + 14, top - 52, bw - 28,
         "A parent and a direction: a maternal loss, a paternal gain.", colour=INK)
    wrap(c, xs[2] + 14, by + 60, bw - 28,
         "One row per chromosome comes out, carrying what it carries, which parent, why in "
         "words, and the counts behind it.")

    for i in (0, 1):
        arrow(c, xs[i] + bw + 6, xs[i + 1] - 6, by + bh / 2)

    # the refusal rail, which is the one thing true of all three steps
    ry, rh = 168, 62
    box(c, xs[0], ry, xs[2] + bw - xs[0], rh, fill=Color(0.96, 0.96, 0.97), edge=RULE)
    text(c, xs[0] + 16, ry + rh - 22, "Any step can stop.", 10.5, INK, "Helvetica-Bold")
    text(c, xs[0] + 125, ry + rh - 22,
         "The row then says which step stopped it and why, in words, instead of going quiet.",
         10.5, INK)
    text(c, xs[0] + 16, ry + 18,
         "That is what a withheld call is: not a missing answer, but a stated reason.", 9, MUTE)
    for i in (0, 1, 2):
        c.setStrokeColor(RULE); c.setLineWidth(0.8); c.setDash(2, 2)
        c.line(xs[i] + bw / 2, by - 4, xs[i] + bw / 2, ry + rh + 4); c.setDash()

    text(c, xs[0], 128, "chr5    one parent's genome    maternal copy absent    "
                        "7% heterozygous over 1,102 markers", 9, INK, "Courier")
    text(c, xs[0], 114, "chr13   no call               withheld                "
                        "genotyped at 37% of its markers", 9, MUTE, "Courier")
    text(c, xs[0], 92, "Two rows of real output, one called and one refused.", 8.5, MUTE)
    c.showPage(); c.save()


out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__),
                                                         "..", "audit", "figures")
os.makedirs(out, exist_ok=True)
fig_dots(os.path.join(out, "share-axis-two-losses.pdf"))
fig_scale(os.path.join(out, "share-scale.pdf"))
fig_schematic(os.path.join(out, "pipeline-three-steps.pdf"))
print("wrote", os.path.abspath(out))
