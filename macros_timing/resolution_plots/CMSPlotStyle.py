import ROOT

""" Global Variables for Styling """
top_margin = 0.07
bottom_margin = 0.12
left_margin = 0.14
right_margin = 0.05

width = 800
height = 800

cmstext = "CMS"
extratext = "Preliminary"

def getStyle():
    """ Get plot style used in CMS
    uses style guidline defined in
    https://twiki.cern.ch/twiki/bin/viewauth/CMS/Internal/FigGuidelines
    especially definitions from
    https://twiki.cern.ch/twiki/pub/CMS/Internal/FigGuidelines/tdrstyle.C
    """

    cmsstyle = ROOT.TStyle('CMS_Style', 'Style for CMS publishing')

    # canvas
    cmsstyle.SetCanvasBorderMode(0)
    cmsstyle.SetCanvasColor(0)
    cmsstyle.SetCanvasDefH(height)
    cmsstyle.SetCanvasDefW(width)
    cmsstyle.SetCanvasDefX(0)
    cmsstyle.SetCanvasDefY(0)

    # pad
    cmsstyle.SetPadBorderMode(0)
    cmsstyle.SetPadColor(0)
    cmsstyle.SetPadGridX(False)
    cmsstyle.SetPadGridY(False)
    cmsstyle.SetGridColor(0)
    cmsstyle.SetGridStyle(3)
    cmsstyle.SetGridWidth(1)

    # frame
    cmsstyle.SetFrameBorderMode(0)
    cmsstyle.SetFrameBorderSize(1)
    cmsstyle.SetFrameFillColor(0)
    cmsstyle.SetFrameFillStyle(0)
    cmsstyle.SetFrameLineColor(1)
    cmsstyle.SetFrameLineStyle(1)
    cmsstyle.SetFrameLineWidth(1)

    # stat box
    cmsstyle.SetOptFile(0)
    cmsstyle.SetOptStat(0)

    # margins
    cmsstyle.SetPadTopMargin(top_margin)
    cmsstyle.SetPadBottomMargin(bottom_margin)
    cmsstyle.SetPadLeftMargin(left_margin)
    cmsstyle.SetPadRightMargin(right_margin)

    # global title
    cmsstyle.SetOptTitle(0)
    cmsstyle.SetTitleFont(42)
    cmsstyle.SetTitleColor(1)
    cmsstyle.SetTitleTextColor(1)
    cmsstyle.SetTitleFillColor(10)
    cmsstyle.SetTitleFontSize(0.05)

    # axis
    cmsstyle.SetAxisColor(1, "XYZ")
    cmsstyle.SetStripDecimals(True)
    cmsstyle.SetTickLength(0.03, "XYZ")
    cmsstyle.SetNdivisions(510, "XYZ")
    cmsstyle.SetPadTickX(1)
    cmsstyle.SetPadTickY(1)

    # axis title
    cmsstyle.SetTitleColor(1, "XYZ")
    cmsstyle.SetTitleFont(42, "XYZ") # Set relative font size
    cmsstyle.SetTitleSize(0.06, "XYZ")
    cmsstyle.SetTitleXOffset(0.9)
    cmsstyle.SetTitleYOffset(1.1)

    # axis labels
    cmsstyle.SetLabelColor(1, "XYZ")
    cmsstyle.SetLabelFont(42, "XYZ") # Set relative font size
    cmsstyle.SetLabelOffset(0.007, "XYZ")
    cmsstyle.SetLabelSize(0.05, "XYZ")

    # postscript
    cmsstyle.SetPaperSize(20.,20.)
    cmsstyle.SetHatchesLineWidth(1)
#    cmsstyle.SetHatchesSpacing(0.05)
    return cmsstyle

def draw_cmstext(position="default",b_pre = True,factor=1):
    x = left_margin + 0.04
    y = 1 - top_margin - top_margin * 0.5
    align =13

    if "default" in position:
        x = left_margin + 0.04
        y = 1 - top_margin - top_margin * 0.5
        align =13
    elif "right" in position:
        x = 1 - right_margin - 0.04
        y = 1 - top_margin - top_margin * 0.5
        align =33

    elif "central" in position:
        x = 1 - right_margin - 0.4
        y = 1 - top_margin - top_margin * 0.5
        align =23

    text1 = ROOT.TLatex(3.5, 24, cmstext)
    text1.SetNDC()
    text1.SetTextAlign(align)
    text1.SetX(x)
    text1.SetY(y)
    text1.SetTextFont(62)
    text1.SetTextSize(0.75*factor *top_margin)
    text1.Draw()


    text2 = ROOT.TLatex(3.5, 24, extratext);
    text2.SetNDC()
    text2.SetTextAlign(align)
    text2.SetX(x)
    text2.SetY(y-top_margin*0.7*factor)
    text2.SetTextFont(52)
    text2.SetTextSize(0.55*text1.GetTextSize())
    if b_pre:
        text2.Draw()

    return [text1,text2]

def draw_info(infotext, x, y, textalign = 31, factor = 1):
    # x = 1 - right_margin - 0.08
    # y = 1 - top_margin + 0.03

    text = ROOT.TLatex(3.5,24, infotext)
    text.SetNDC()
    text.SetTextAlign(textalign)
    text.SetX(x)
    text.SetY(y)
    text.SetTextFont(42)
    text.SetTextSize(0.6*top_margin*factor)
    text.Draw()

    return text

def draw_lumi(simonly =False, offset = 0.2, lumi="35.9"):
    cmsTextSize      = 0.75
    extraOverCmsTextSize  = 0.76
    lumiTextSize     = 0.7
    lumiTextOffset   = offset
    cmsTextFont   = 61

    lumiText=lumi+" fb^{-1} (13 TeV)"
    if "nb" in lumi: lumiText = lumi + " (13 TeV)"
    if simonly:    lumiText="(13 TeV)"

    latex = ROOT.TLatex()
    latex.SetNDC();
    latex.SetTextAngle(0);
    latex.SetTextColor(ROOT.kBlack);

    extraTextSize = extraOverCmsTextSize*cmsTextSize;

    latex.SetTextFont(42);
    latex.SetTextAlign(31);
    latex.SetTextSize(lumiTextSize*top_margin);
    latex.DrawLatex(1-right_margin,1-top_margin+lumiTextOffset*top_margin,lumiText);




    return latex
