import ROOT
ROOT.gROOT.SetBatch()

primaryKaons   = lambda x: 'motherindex != -1 && abs(catchall) == '+str(x)
secondaryKaons = lambda x: 'motherindex == -1 && abs(catchall) == '+str(x)
Kaons = lambda x: 'abs(catchall) == '+str(x)+' && (abs(pdgId) != 3 || abs(pdgId) != 2)'

f = ROOT.TFile('hzkin.root')
t = f.Get('mctrue')

plots = [ ('p',''),('cos(theta)',''),('p:cos(theta)','COLZ') ]

for res in [ 23, 25 ]:
    c = ROOT.TCanvas()
    c.Divide(3,4)

    for i,(var,opt) in enumerate(plots):
        c.cd(3*i+1)
        t.Draw(var,primaryKaons(res),opt)
        
        c.cd(3*i+2)
        t.Draw(var,secondaryKaons(res),opt)
    
        c.cd(3*i+3)
        t.Draw(var,Kaons(res),opt)
    
    c.cd(3*i+4)
    
    t.Draw('p>>hpri(100,0,60)',primaryKaons(res))
    hpri = ROOT.gDirectory.Get('hpri')
    hpri.Scale(1./hpri.Integral())
    #hpri.Sumw2()

    t.Draw('p>>hsec(100,0,60)',secondaryKaons(res))
    hsec = ROOT.gDirectory.Get('hsec')
    hsec.Scale(1./hsec.Integral())
    #hsec.Sumw2()

    h = hpri.Clone('dividedh')
    #h.Sumw2()
    h.Divide(hsec)

    h.Draw()#"PE")


    c.SaveAs('kaonsat'+str(res)+'.pdf')




