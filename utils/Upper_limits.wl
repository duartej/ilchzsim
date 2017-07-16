(* ::Package:: *)

Quit[]


(* ::Section:: *)
(*Read efficiencies and calculate the upper limits on the signal strength*)


(* ::Subsection:: *)
(*Input*)


experiment="CLIC";
sqrts="350";


If[experiment=="CLIC",
	Which[
	sqrts=="350",
	{
		Lumi=500.;
		process="HZ";
		NHiggs["inv"]=(8000./brbbbar + 372./brccbar)/2;
		NHiggs["had"]=(11100./brbbbar+ 434./brccbar)/2;
		Nbkg["inv"]=2100+2090+104+30+1230;
		Nbkg["had"]=60+89+9990+11400;
		Zmodes={"inv","had"};
	},
	sqrts=="1400",
	{
		Lumi=1500.;
		process="Hinv";
		NHiggs["none"]=(65400./brbbbar + 3790./brccbar)/2;
		Nbkg["none"]=18500+23600+18500+170000+22200;
		Zmodes={"none"};
	},
	sqrts=="3000",
	{
		Lumi=2000.;
		process="Hinv";
		NHiggs["none"]=(120000./brbbbar + 6380./brccbar)/2;
		Nbkg["none"]=47400+52200+118000+394000+207000;
		Zmodes={"none"};
	}
	]
,
Print["did not yet implement experiment "<>experiment]
];


channel="CC";


(* ::Subsubsection::Closed:: *)
(*Input directory*)


SetDirectory["/home/matthias/Physics/Projects/ssbar/code/runcode/2sqrt2_"<>sqrts<>"_"<>process];


(* ::Subsubsection::Closed:: *)
(*Branching ratios*)


brbbbar =5.66 10^-1;
brccbar = 2.85 10^-2;
brssbar = 2.41 10^-4;
brgg = 8.50 10^-2;
bruubar = 1.3 10^-7;
brddbar = 5.8 10^-8;

flavors={"b","c","s","u","d","g"};
bkgflavors={"b","c","u","d","g"};
BR[mode_]:=Which[
	mode=="b",brbbbar,
	mode=="c",brccbar,
	mode=="s",brssbar,
	mode=="u",bruubar,
	mode=="d",brddbar,
	mode=="g",brgg];


(* ::Subsubsection::Closed:: *)
(*Load Colors*)


col1=ColorData[81,"ColorList"]
col2=ColorData[112,"ColorList"]


(* ::Subsection:: *)
(*Functions*)


(* ::Subsubsection::Closed:: *)
(*functions for info extraction and event calculation*)


GetInfos[filelist_]:=Module[{},
	Return[Map[StringSplit[#,{"_","-"}][[2;;6]]&,filelist]]
	];


GetEfficiencyFunction[flavor_,efftable_]:=Module[{column,f},
	column=Which[
		flavor=="b",2,
		flavor=="c",3,
		flavor=="s",4,
		flavor=="u",5,
		flavor=="d",6,
		flavor=="g",7];
	f=Interpolation[Map[{#[[1]],#[[column]]}&,efftable]];
	Return[f];
];	


eff[d0_,effK_,flavor_]:=GetEfficiencyFunction[flavor,efficiencytable[d0,effK]];


evts[d0_,effK_,flavor_,Zmode_,pcut_]:=NHiggs[Zmode]*BR[flavor]*eff[d0,effK,flavor][pcut]


(* signal and background (only Higgs) numbers *)
SBHiggs[d0_,effK_,Zmode_,pcut_]:={evts[d0,effK,"s",Zmode,pcut],Map[evts[d0,effK,#,Zmode,pcut]&,bkgflavors]//Total}
significanceHiggs[d0_,effK_,Zmode_,pcut_]:=Map[#[[1]]/Sqrt[#[[2]]]&,{SBHiggs[d0,effK,Zmode,pcut]}][[1]]


(* signal and all background numbers *)
SBnoHiggs[d0_,effK_,Zmode_,pcut_]:={0,Nbkg[Zmode]*(eff[d0,effK,"u"][pcut]+eff[d0,effK,"d"][pcut]+eff[d0,effK,"s"][pcut]+eff[d0,effK,"c"][pcut]+eff[d0,effK,"b"][pcut])/5};
SBall[d0_,effK_,Zmode_,pcut_]:=SBHiggs[d0,effK,Zmode,pcut]+SBnoHiggs[d0,effK,Zmode,pcut];
significanceall[d0_,effK_,Zmode_,pcut_]:=Map[#[[1]]/Sqrt[#[[2]]]&,{SBall[d0,effK,Zmode,pcut]}][[1]]


(* ::Subsubsection:: *)
(*Limit related functions*)


(* ::Text:: *)
(*Based on PDG statistics review 2016*)


sup[b_,nobs_,CLgoal_]:=1/2 Quantile[ChiSquareDistribution[2(nobs+1)],1-(1-CLgoal)*(1-CDF[ChiSquareDistribution[2(nobs+1)],2*b])]-b


(* Upper limit on the signal strength at the 
	CLgoal confidence level
when there are 
	s signal events
	b background events
	nobs events observed.
*)

muUp[s_,b_,CLgoal_]:=Module[{},
	Return[sup[b,s+b,CLgoal]/s]
];


(* ::Subsection:: *)
(*Read file and get upper limits*)


efffiles=FileNames["*txt",{channel}];
efficiencies=GroupBy[GetInfos[efffiles],First];


d0cut=Keys[efficiencies];

For[i=1, i<=Length[d0cut], i++,
	d0=d0cut[[i]];
	pids=efficiencies[d0];
	PrintTemporary["analyzing d0cut="<>d0];
	
	d0Efficiencies[d0]={};
	For[j=1, j<=Length[pids], j++,
		pid=efficiencies[d0][[j]];
		effK=pid[[3]];
		effP=pid[[4]];
		AppendTo[d0Efficiencies[d0],effK];
		
		efficiencytable[d0, effK]=Select[Import[channel<>"/efficiencies_"<>d0<>"_0.95-"<>effK<>"-"<>effP<>"-0.9-PID.txt","Table"], Length[#]>0&];
		For[k=1,k<=Length[Zmodes], k++,
			Zmode=Zmodes[[k]];
			ULH[d0,effK,Zmode]=Interpolation[Map[{#, muUp[SBHiggs[d0,effK,Zmode,#][[1]],SBHiggs[d0,effK,Zmode,#][[2]],0.95]}&,Range[10,20,1]]];
			ULA[d0,effK,Zmode]=Interpolation[Map[{#, muUp[SBall[d0,effK,Zmode,#][[1]],SBall[d0,effK,Zmode,#][[2]],0.95]}&,Range[10,20,1]]];
			UpperLimitHiggs[d0,effK,Zmode,pcut_]=ULH[d0,effK,Zmode][pcut];
			UpperLimitall[d0,effK,Zmode,pcut_]=ULA[d0,effK,Zmode][pcut];
		];
	];
];


SBall["0.025","1.0","had",10]


sup[3182,3183,0.95]


(* ::Subsection:: *)
(*Plot*)


For[i=1, i<=Length[d0cut], i++,
	d0=d0cut[[i]];
	Print["\!\(\*SubscriptBox[\(d\), \(0\)]\) cut = "<>d0];
	For[j=1, j<=Length[d0Efficiencies[d0]], j++,
		effK=d0Efficiencies[d0][[j]];
	
		effPlot[d0,effK]=LogPlot[{eff[d0,effK,"s"][x],
							eff[d0,effK,"b"][x],
							eff[d0,effK,"c"][x],
							eff[d0,effK,"u"][x],
							eff[d0,effK,"d"][x],
							eff[d0,effK,"g"][x]},{x,10,20},
			Frame->True,
			FrameLabel->{"\!\(\*SubsuperscriptBox[\(p\), \(||\), \(cut\)]\) [GeV]","\[Epsilon]"},
			BaseStyle->13,
			AspectRatio->1,
			ImageSize->300,
			PlotLabel->"\!\(\*SubscriptBox[\(d\), \(0\)]\)="<>d0<>"mm, \!\(\*SubscriptBox[\(\[Epsilon]\), \(K\)]\)="<>effK,
			PlotLegends->{"s","b","c","u","d","g"}];

		For[k=1, k<=Length[Zmodes], k++,
			Zmode=Zmodes[[k]];
			eventsPlot[d0,effK,Zmode]=LogPlot[{evts[d0,effK,"s",Zmode,x],
									evts[d0,effK,"b",Zmode,x],
									evts[d0,effK,"c",Zmode,x],
									evts[d0,effK,"u",Zmode,x],
									evts[d0,effK,"d",Zmode,x],
									evts[d0,effK,"g",Zmode,x](*,
									SBnoHiggs[d0,effK,Zmode,x][[2]]*)},{x,10,20},
				Frame->True,
				FrameLabel->{"\!\(\*SubsuperscriptBox[\(p\), \(||\), \(cut\)]\) [GeV]","# events per channel"},
				BaseStyle->13,
				AspectRatio->1,
				ImageSize->300,
				PlotLabel->"\!\(\*SubscriptBox[\(d\), \(0\)]\)="<>d0<>"mm, \!\(\*SubscriptBox[\(\[Epsilon]\), \(K\)]\)="<>effK<>", Z\[Rule]"<>Zmode,
				PlotLegends->{"s","b","c","u","d","g","non-Higgs"}];
		];

		If[d0=="0.017",
			p=Row[{effPlot[d0,effK], eventsPlot[d0,effK,Zmodes[[1]]],eventsPlot[d0,effK,Zmodes[[2]]]},"         "];
			Print[p];
		];
	];
];


sigRange=Full;
Which[
	sqrts=="350",{sigRange={0,0.2}}
];

For[i=1, i<=Length[d0cut], i++,
	d0=d0cut[[i]];
	Print["\!\(\*SubscriptBox[\(d\), \(0\)]\) cut = "<>d0];

	For[j=1,j<=Length[Zmodes],j++,
		legend=LineLegend[{col1[[1]],col1[[2]],col1[[3]]},Map["\!\(\*SubscriptBox[\(\[Epsilon]\), \(K\)]\)="<>#&,d0Efficiencies[d0]],LegendLayout->"Column"];
		significancePlotHiggs[d0,Zmodes[[j]]]=Plot[Map[significanceHiggs[d0,#,Zmodes[[j]],x]&,d0Efficiencies[d0]],{x,10,20},
			PlotRange->{Full,sigRange},
			Frame->True,
			FrameLabel->{"\!\(\*SubsuperscriptBox[\(p\), \(||\), \(cut\)]\) [GeV]","S/\!\(\*SqrtBox[\(B\)]\)"},
			PlotLabel->"\!\(\*SubsuperscriptBox[\(d\), \(0\), \(cut\)]\)="<>d0<>"mm,  Z \[Rule] "<>Zmodes[[j]]<>", Higgs only",
			PlotStyle->{col1[[1]],col1[[2]],col1[[3]]},
			BaseStyle->13,
			ImageSize->300,
			Evaluated->True,
			PlotLegends->Placed[legend,{0.1,0.6}]];
		significancePlotAll[d0,Zmodes[[j]]]=Plot[Map[significanceall[d0,#,Zmodes[[j]],x]&,d0Efficiencies[d0]],{x,10,20},
			PlotRange->{Full,sigRange},
			Frame->True,
			FrameLabel->{"\!\(\*SubsuperscriptBox[\(p\), \(||\), \(cut\)]\) [GeV]","S/\!\(\*SqrtBox[\(B\)]\)"},
			PlotLabel->"\!\(\*SubsuperscriptBox[\(d\), \(0\), \(cut\)]\)="<>d0<>"mm,  Z \[Rule] "<>Zmodes[[j]],
			PlotStyle->{col1[[1]],col1[[2]],col1[[3]]},
			BaseStyle->13,
			ImageSize->300,
			Evaluated->True,
			PlotLegends->Placed[legend,{0.5,0.8}]];
		UpperLimitPlotHiggs[d0,Zmodes[[j]]]=Plot[Map[UpperLimitHiggs[d0,#,Zmodes[[j]],x]&,d0Efficiencies[d0]],{x,10,20},
			Frame->True,
			FrameLabel->{"\!\(\*SubsuperscriptBox[\(p\), \(||\), \(cut\)]\) [GeV]","95% CL on \[Mu]"},
			PlotLabel->"\!\(\*SubsuperscriptBox[\(d\), \(0\), \(cut\)]\)="<>d0<>"mm,  Z \[Rule] "<>Zmodes[[j]]<>", Higgs only",
			PlotStyle->{col1[[1]],col1[[2]],col1[[3]]},
			BaseStyle->13,
			ImageSize->300,
			Evaluated->True,
			AspectRatio->1,
			PlotLegends->Placed[legend,{0.25,0.7}]];
		UpperLimitPlotAll[d0,Zmodes[[j]]]=Plot[Map[UpperLimitall[d0,#,Zmodes[[j]],x]&,d0Efficiencies[d0]],{x,10,20},
			Frame->True,
			FrameLabel->{"\!\(\*SubsuperscriptBox[\(p\), \(||\), \(cut\)]\) [GeV]","95% CL on \[Mu]"},
			PlotLabel->"\!\(\*SubsuperscriptBox[\(d\), \(0\), \(cut\)]\)="<>d0<>"mm,  Z \[Rule] "<>Zmodes[[j]],
			PlotStyle->{col1[[1]],col1[[2]],col1[[3]]},
			BaseStyle->13,
			ImageSize->300,
			Evaluated->True,
			AspectRatio->1,
			PlotLegends->Placed[legend,{0.1,0.7}]];
		Print[Grid[{{significancePlotHiggs[d0,Zmodes[[j]]],significancePlotAll[d0,Zmodes[[j]]]},{UpperLimitPlotHiggs[d0,Zmodes[[j]]],UpperLimitPlotAll[d0,Zmodes[[j]]]}}]];
	];
];


pcut=10;
For[k=1, k<=Length[Zmodes], k++,
	Zmode=Zmodes[[k]];
	Print[ListLinePlot[
		{Map[{ToExpression[#],UpperLimitall[#,"0.8",Zmode,10]}&,d0cut],
		Map[{ToExpression[#],UpperLimitall[#,"0.95",Zmode,10]}&,d0cut],
		Map[{ToExpression[#],UpperLimitall[#,"1.0",Zmode,10]}&,d0cut]},
		Frame->True,
		FrameLabel->{"\!\(\*SubsuperscriptBox[\(d\), \(0\), \(cut\)]\) [mm]", "95% CL on \[Mu]"},
		PlotRange->{Full,{50,Automatic}},
		PlotLegends->{"0.8","0.95","1.0"},
		BaseStyle->13,
		PlotLabel->"Z \[Rule] "<>Zmode<>", \!\(\*SubsuperscriptBox[\(p\), \(||\), \(cut\)]\) = "<>ToString[pcut]<>" GeV"]];
	];



