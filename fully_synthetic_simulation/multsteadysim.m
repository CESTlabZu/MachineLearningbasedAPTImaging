%p3 = runsteadysim1(6.983, .122, 1, r1s, r2s, r1w, r2w, 0, .03, gauss, 500, .05, 140, 0, .008, .002, 7)
%p4 = runsteadysim1(6.983, .122, 1, r1s, r2s, r1w, r2w, 0, .03, gauss, 1000, .05, 140, 0, .008, .002, 7)
%p5 = runsteadysim1(6.983, .122, 1, r1s, r2s, r1w, r2w, 0, .03, gauss, 500, .05, 140, 1, .008, .002, 7)
%p6 = runsteadysim1(6.983, .122, 1, r1s, r2s, r1w, r2w, 0, .03, gauss, 1000, .05, 140, 1, .008, .002, 7)
%p7 = runsteadysim1(6.983, .122, 1, r1s, r2s, r1w, r2w, 0, .01, gauss, 300, .025, 140, 0, .008, .002, 7)
%p8 = runsteadysim1(6.983, .122, 1, r1s, r2s, r1w, r2w, 0, .01, gauss, 500, .025, 140, 0, .008, .002, 7)
%p9 = runsteadysim1(6.983, .122, 1, r1s, r2s, r1w, r2w, 0, .01, gauss, 300, .025, 140, 1, .008, .002, 7)
%p10 = runsteadysim1(6.983, .122, 1, r1s, r2s, r1w, r2w, 0, .01, gauss, 500, .025, 140, 1, .008, .002, 7)
%p11 = runsteadysim1(6.983, .122, 1, r1s, r2s, r1w, r2w, 0, .01, gauss, 700, .025, 140, 0, .008, .002, 7)
%p12 = runsteadysim1(6.983, .122, 1, r1s, r2s, r1w, r2w, 0, .01, gauss, 700, .025, 140, 1, .008, .002, 7)
%p13 = runsteadysim1(6.983, .122, 1, r1s, r2s, r1w, r2w, 0, .03, gauss, 1400, .05, 140, 0, .008, .002, 7)
%p14 = runsteadysim1(6.983, .122, 1, r1s, r2s, r1w, r2w, 0, .03, gauss, 1400, .05, 140, 1, .008, .002, 7)
%p7b = runsteadysim2(6.983, .122, 1, r1s, r2s, r1w, r2w, 0, .01, gauss, 300, .025, 140, 0, .008, .002, 7)
%p9b = runsteadysim2(6.983, .122, 1, r1s, r2s, r1w, r2w, 0, .01, gauss, 300, .025, 140, 1, .008, .002, 7)
%p8b = runsteadysim2(6.983, .122, 1, r1s, r2s, r1w, r2w, 0, .01, gauss, 500, .025, 140, 0, .008, .002, 7)
%p10b = runsteadysim2(6.983, .122, 1, r1s, r2s, r1w, r2w, 0, .01, gauss, 500, .025, 140, 1, .008, .002, 7)


%suncest1 = runsteadysim1(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .02, gauss, 2*ut, .05, 140, 0, .002, .001, 7)
%suncest2 = runsteadysim1(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .025, gauss, 2*ut, .05, 140, 0, .002, .001, 7)
%suncest3 = runsteadysim1(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .03, gauss, 2*ut, .05, 140, 0, .002, .001, 7)
%suncest4 = runsteadysim1(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 2*ut, .05, 140, 0, .002, .001, 7)
%suncest5 = runsteadysim1(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .04, gauss, 2*ut, .05, 140, 0, .008, .002, 7)
%suncestexc1 = runsteadysim1(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .02, gauss, 2*ut, .05, 140, 1, .002, .001, 7);
%suncestexc2 = runsteadysim1(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .025, gauss, 2*ut, .05, 140, 1, .002, .001, 7);
%suncestexc3 = runsteadysim1(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .03, gauss, 2*ut, .05, 140, 1, .002, .001, 7);
%suncestexc4 = runsteadysim1(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 2*ut, .05, 140, 1, .002, .001, 7);
%suncestexc5 = runsteadysim1(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .04, gauss, 2*ut, .05, 140, 1, .007, .002, 7);

%suncest1p = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 1.25*ut, .05, 140, 0, .002, .001, 7,l)
%suncest2p = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 1.3*ut, .05, 140, 0, .002, .001, 7,l)
%suncest3p = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 1.35*ut, .05, 140, 0, .002, .001, 7,l)
%suncest4p = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 1.4*ut, .05, 140, 0, .002, .001, 7,l)
%suncest5p = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 1.45*ut, .05, 140, 0, .002, .001, 7, l)
%suncest6p = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 1.5*ut, .05, 140, 0, .002, .001, 7,l)
%suncest7p = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 1.55*ut, .05, 140, 0, .002, .001, 7,l)
%suncest8p = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 1.6*ut, .05, 140, 0, .002, .001, 7,l)
%suncest9p = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 1.65*ut, .05, 140, 0, .002, .001, 7,ll)
%suncest10p = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 1.7*ut, .05, 140, 0, .002, .001, 7,ll)
%suncest11p = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 1.9*ut, .05, 140, 0, .002, .001, 7,ll)
%suncest12p = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 2.1*ut, .05, 140, 0, .002, .001, 7,ll)
%suncest13p = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 3*ut, .05, 140, 0, .002, .001, 7,ll)
%suncest14p = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 3.5*ut, .05, 140, 0, .002, .001, 7,ll)
%suncest15p = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 5*ut, .05, 140, 0, .002, .001, 7,l)

%suncest3p = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 2.5*ut, .05, 140, 0, .002, .001, 7,ll)
%suncest4p = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 2.75*ut, .05, 140, 0, .002, .001, 7,ll)
%suncest5p = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 3*ut, .05, 140, 0, .002, .001, 7,ll)
%suncest6p = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 3.25*ut, .05, 140, 0, .002, .001, 7,ll)
%suncest7p = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 3.5*ut, .05, 140, 0, .002, .001, 7,ll)
%suncest8p = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 4*ut, .05, 140, 0, .002, .001, 7,ll)
%suncest9p = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 4.5*ut, .05, 140, 0, .002, .001, 7,ll)
%suncest11p = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 5.5*ut, .05, 140, 0, .002, .001, 7,ll)
%suncest12p = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 6*ut, .05, 140, 0, .002, .001, 7,ll)
%suncest1p = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 1.5*ut, .05, 140, 0, .002, .001, 7,ll)
%suncest2pa = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 1.7893*ut, .05, 140, 0, .002, .001, 7,ll)
%suncest2pb = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .03, gauss, 1.7893*ut, .05, 140, 0, .002, .001, 7,ll)
%suncest2pc1 = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .015, gauss, 1.5*ut, .05, 140, 0, .002, .001, 7,ll)
%suncest2pc1a = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .012, gauss, 1.5*ut, .05, 140, 0, .002, .001, 7,ll)
%suncest2pc1b = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .03, gauss, 1.5*ut, .05, 140, 0, .002, .001, 7,ll)

%suncestpulsedw1t1 = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 1.5*ut, .05, 140, 0, .002, .001, 7,ll);
%suncestpulsedw1t2 = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 1.75*ut, .05, 140, 0, .002, .001, 7,ll);
%suncestpulsedw1t3 = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 2.0*ut, .05, 140, 0, .002, .001, 7,ll);
%suncestpulsedw1t4 = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 2.25*ut, .05, 140, 0, .002, .001, 7,ll);
%suncestpulsedw1t5 = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 2.5*ut, .05, 140, 0, .002, .001, 7,ll);

%suncestpulsed1s = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .015, gauss, 1.5*ut, .05, 140, 0, .002, .001, 7,ll)
%suncestpulsed2s = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .025, gauss, 1.5*ut, .05, 140, 0, .002, .001, 7,ll)
%suncestpulsed3s = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 1.5*ut, .05, 140, 0, .002, .001, 7,ll)
%suncestpulsed3t = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 1.78*ut, .05, 140, 0, .002, .001, 7,pp)
%suncestpulsed4s = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .045, gauss, 1.5*ut, .05, 140, 0, .002, .001, 7,ll)

%suncestpulsed4s1 = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .049, gauss, 1.5*ut, .05, 140, 0, .002, .001, 7,ll)

%suncestpulsedw1 = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .015, gauss, 1.5*ut, .05, 140, 0, .002, .001, 7,ll)
%suncestpulsedw2 = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .015, gauss, 1.75*ut, .05, 140, 0, .002, .001, 7,ll)
%suncestpulsedw3 = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .015, gauss, 2.0*ut, .05, 140, 0, .002, .001, 7,ll)
%suncestpulsedw4 = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .015, gauss, 2.25*ut, .05, 140, 0, .002, .001, 7,ll)

%suncestasym1s = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .025, gauss, 1.5*ut, .1, 140, 0, .002, .001, 7,ll)
%suncestasym2s = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .025, gauss, 1.5*ut, .05, 140, 0, .002, .001, 7,ll)
%suncestasym3s = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .025, gauss, 1.5*ut, .15, 140, 0, .002, .001, 7,ll)
%suncestasym4s = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .025, gauss, 1.5*ut, .06, 140, 0, .002, .001, 7,ll)
%suncestasym5s = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .025, gauss, 1.5*ut, .08, 140, 0, .002, .001, 7,ll)

%suncest2pc2 = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .025, gauss, 1.5*ut, .05, 140, 0, .002, .001, 7,ll)
%suncest2pc3 = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 1.5*ut, .05, 140, 0, .002, .001, 7,ll)
%suncest2pc4 = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .045, gauss, 1.5*ut, .05, 140, 0, .002, .001, 7,ll)
%suncest3p = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 2.0*ut, .05, 140, 0, .002, .001, 7,ll)
%suncest5p1 = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 1.455*ut, .05, 140, 0, .002, .001, 7, l);
%suncest5p2 = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 1.46*ut, .05, 140, 0, .002, .001, 7, l);
%suncest5p3 = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 1.465*ut, .05, 140, 0, .002, .001, 7, l);
%suncest5p4 = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 1.47*ut, .05, 140, 0, .002, .001, 7, l);

%suncestexc1c = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .02, gauss, ut, .05, 140, 1, .002, .001, 7,k);
%suncestexc2a = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .025, gauss, 2*ut, .05, 140, 1, .002, .001, 7,k);
%suncestexc3a = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .03, gauss, 2*ut, .05, 140, 1, .002, .001, 7,k);
%suncestexc4c = runsteadysim2(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, ut, .05, 140, 1, .002, .001, 7,k);


%sunparacest1 = runsteadysim1(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .02, gauss, 2*ut, .05, 140, 0, .008, .002, 7)
%sunparacest2 = runsteadysim1(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .025, gauss, 2*ut, .05, 140, 0, .008, .002, 7)
%sunparacest3 = runsteadysim1(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .03, gauss, 2*ut, .05, 140, 0, .008, .002, 7)
%sunparacest4 = runsteadysim1(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 2*ut, .05, 140, 0, .008, .002, 7, k)
%sunparacest5 = runsteadysim1(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .04, gauss, 2*ut, .05, 140, 0, .008, .002, 7, k)
%sunparacestexc1 = runsteadysim1(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .02, gauss, 2*ut, .05, 140, 1, .007, .002, 7);
%sunparacestexc2 = runsteadysim1(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .025, gauss, 2*ut, .05, 140, 1, .007, .002, 7);
%sunparacestexc3 = runsteadysim1(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .03, gauss, 2*ut, .05, 140, 1, .007, .002, 7);
%sunparacestexc4 = runsteadysim1(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 2*ut, .05, 140, 1, .007, .002, 7, k);
%sunparacestexc5 = runsteadysim1(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .04, gauss, 2*ut, .05, 140, 1, .007, .002, 7, k);

%gausspl1 = runsteadysimgauss(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .025, gauss, 1.5*ut, .05, 140, 0, .002, .001, 7, k, 0);
%gausspl2 = runsteadysimgauss(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .030, gauss, 1.5*ut, .05, 140, 0, .002, .001, 7, k, 0);
gausspl3 = runsteadysimgauss(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 1.5*ut, .05, 140, 0, .002, .001, 7, k, 0);

%gaussexcpl1 = runsteadysimgauss(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .025, gauss, 1.5*ut, .05, 140, 1, .002, .001, 7, k, 0);
%gaussexcpl2 = runsteadysimgauss(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .030, gauss, 1.5*ut, .05, 140, 1, .002, .001, 7, k, 0);
%gaussexcpl3 = runsteadysimgauss(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 1.5*ut, .05, 140, 1, .002, .001, 7, k, 0);

squarepl1 = runsteadysimsquare(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .025, gauss, 1.5*ut, .05, 140, 0, .002, .001, 7, k, 0);
%squarepl2 = runsteadysimsquare(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .030, gauss, 1.5*ut, .05, 140, 0, .002, .001, 7, k, 0);
%squarepl3 = runsteadysimsquare(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 1.5*ut, .05, 140, 0, .002, .001, 7, k, 0);

%squareexcpl1 = runsteadysimsquare(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .025, gauss, 1.5*ut, .05, 140, 1, .002, .001, 7, k, 0);
%squareexcpl2 = runsteadysimsquare(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .030, gauss, 1.5*ut, .05, 140, 1, .002, .001, 7, k, 0);
%squareexcpl3 = runsteadysimsquare(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 1.5*ut, .05, 140, 1, .002, .001, 7, k, 0);

%gausspl1f = runsteadysimgauss(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .025, gauss, 1.5*ut, .05, 140, 0, .002, .001, 7, k, 1);
%gausspl2s = runsteadysimgauss(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .030, gauss, 1.5*ut, .05, 140, 0, .002, .001, 7, k, 1);
%gausspl3s = runsteadysimgauss(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 1.5*ut, .05, 140, 0, .002, .001, 7, k, 1);

%gaussexcpl1s = runsteadysimgauss(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .025, gauss, 1.5*ut, .05, 140, 1, .002, .001, 7, k, 1);
%gaussexcpl2s = runsteadysimgauss(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .030, gauss, 1.5*ut, .05, 140, 1, .002, .001, 7, k, 1);
%gaussexcpl3s = runsteadysimgauss(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 1.5*ut, .05, 140, 1, .002, .001, 7, k, 1);

%squarepl1s = runsteadysimsquare(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .025, gauss, 1.5*ut, .05, 140, 0, .002, .001, 7, k, 1);
%squarepl2s = runsteadysimsquare(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .030, gauss, 1.5*ut, .05, 140, 0, .002, .001, 7, k, 1);
%squarepl3s = runsteadysimsquare(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 1.5*ut, .05, 140, 0, .002, .001, 7, k, 1);

%squareexcpl1s = runsteadysimsquare(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .025, gauss, 1.5*ut, .05, 140, 1, .002, .001, 7, k, 1);
%squareexcpl2s = runsteadysimsquare(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .030, gauss, 1.5*ut, .05, 140, 1, .002, .001, 7, k, 1);
%squareexcpl3s = runsteadysimsquare(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 1.5*ut, .05, 140, 1, .002, .001, 7, k, 1);

%gausspl1s = runsteadysimgauss(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .0066, gauss, 50, .01, 140, 0, .002, .001, 7, k, 1);
%gausspl2s = runsteadysimgauss(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .030, gauss, 1.5*ut, .05, 140, 0, .002, .001, 7, k, 1);
%gausspl3s = runsteadysimgauss(50, .005, 1, 1, 1/.015, 1/3, 1/.05, 760, .035, gauss, 1.5*ut, .05, 140, 0, .002, .001, 7, k, 1);
