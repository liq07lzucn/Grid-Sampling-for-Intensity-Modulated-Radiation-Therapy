function [caseParams oMap] = getCaseParams(patID)
%function to encapsulate case-specific parameter definition
global stateS
eval(['init_pat_', int2str(patID)]);

caseParams.patID = patID;


try
    caseParams.nameRun = nameRun;
    
    organsMap = containers.Map();
    if patID == 1
        caseParams.PTV1 = PTV1;
        caseParams.PTV1wMoat = PTV1wMoat; % 3D expansion of the target, say by 0.5 cm.
        caseParams.PTV2 = PTV2;
        caseParams.PTV2wMoat = PTV2wMoat;
        caseParams.skin = skin;
        
        caseParams.Brain = Brain;
        caseParams.BrainStem = BrainStem;
        caseParams.Mandible = Mandible;
        
        caseParams.LeftOrbit = LeftOrbit;
        caseParams.OpticChiasm = OpticChiasm;
        caseParams.OpticNerves = OpticNerves;
        
        caseParams.PatrodGland = PatrodGland;
        caseParams.RightOrbit = RightOrbit;
        caseParams.SpinalCord = SpinalCord;
        caseParams.PTV3 = PTV3;
        caseParams.PTV3wMoat = PTV3wMoat;
        
        organsMap(int2str(PTV1)) = 'Target1';
        organsMap(int2str(PTV2)) = 'Target2';
        organsMap(int2str(PTV3)) = 'PTV3';
        organsMap(int2str(skin)) = 'skin';
        
        
        organsMap(int2str(Brain)) = 'Brain';
        organsMap(int2str(BrainStem)) = 'BrainStem';
        organsMap(int2str(Mandible)) = 'Mandible';
        organsMap(int2str(LeftOrbit)) = 'LeftOrbit';
        organsMap(int2str(OpticNerves)) = 'OpticNerves';
        
        organsMap(int2str(PatrodGland)) = 'PatrodGland';
        organsMap(int2str(OpticChiasm)) = 'OpticChiasm';
        
        organsMap(int2str(RightOrbit)) = 'RightOrbit';
        organsMap(int2str(SpinalCord)) = 'SpinalCord';
        
        
        oMap = organsMap;
    elseif patID==2
        caseParams.PTV1 = PTV1;
        caseParams.PTV1wMoat = PTV1wMoat; % 3D expansion of the target, say by 0.5 cm.
        caseParams.PTV2 = PTV2;
        caseParams.PTV2wMoat = PTV2wMoat;
        caseParams.skin = skin;
        caseParams.Rectum = RECTUM;
        caseParams.Bladder = Bladder;
        caseParams.FemurLeft = FEMUUR_LEFT;
        
        caseParams.FemurRight = FEMUUR_RIGHT;
        caseParams.PTV3 = GTV;
        
        organsMap(int2str(PTV1)) = 'Target1';
        organsMap(int2str(PTV2)) = 'Target2';
        organsMap(int2str(PTV3)) = 'PTV3';
        organsMap(int2str(skin)) = 'skin';
        
        
        organsMap(int2str(Rectum)) = 'Rectum';
        organsMap(int2str(Bladder)) = 'Bladder';
        organsMap(int2str(FemurLeft)) = 'FemurLeft';
        organsMap(int2str(FemurRight)) = 'FemurRight';
        organsMap(int2str(GTV)) = 'GTV';
        oMap = organsMap;
    elseif patID ==3
        caseParams.PTV1 = PTV1;
        caseParams.PTV1wMoat = PTV1wMoat; % 3D expansion of the target, say by 0.5 cm.
        caseParams.PTV2 = PTV2;
        caseParams.PTV2wMoat = PTV2wMoat;
        caseParams.SKIN = SKIN;
        caseParams.CTV1 = CTV1;
        caseParams.ESOPHAGUS = ESOPHAGUS;
        caseParams.HEART = HEART;
        caseParams.LIVER = LIVER;
        caseParams.LUNG_CONTRA = LUNG_CONTRA;
        caseParams.LUNG_IPSI = LUNG_IPSI;
        caseParams.SPINAL_CORD = SPINAL_CORD;
        caseParams.INITIAL_REF = INITIAL_REF;
        caseParams.TOTAL_LUNG = TOTAL_LUNG;
        caseParams.PTV3 = GTV1;
        
        organsMap(int2str(PTV1)) = 'Target1';
        organsMap(int2str(PTV2)) = 'Target2';
        organsMap(int2str(GTV1)) = 'PTV3';
        organsMap(int2str(SKIN)) = 'skin';
        organsMap(int2str(SPINAL_CORD)) = 'SPINAL_CORD';
        organsMap(int2str(LUNG_IPSI)) = 'LUNG_IPSI';
        organsMap(int2str(LUNG_CONTRA)) = 'LUNG_CONTRA';
        organsMap(int2str(LIVER)) = 'LIVER';
        organsMap(int2str(HEART)) = 'HEART';
        organsMap(int2str(ESOPHAGUS)) = 'ESOPHAGUS';
        organsMap(int2str(CTV1)) = 'CTV1';
        organsMap(int2str(INITIAL_REF)) = 'INITIAL_REF';
        organsMap(int2str(TOTAL_LUNG)) = 'TOTAL_LUNG';
        oMap = organsMap;
    elseif patID == 5
        caseParams.Skin = Skin;
        caseParams.PTV1 = PTV1;
        caseParams.Brain = Brain;
        caseParams.BrainStem = BrainStem;
        caseParams.OpticNerves = OpticNerves;
        caseParams.Eyes = Eyes;
        caseParams.Parotids = Parotids;
        caseParams.Cord = Cord;
        caseParams.PTV2 = PTV2;
        caseParams.N2A = N2A;
        caseParams.PTV3 = PTV3;
        
        organsMap(int2str(Skin)) = 'skin';
        organsMap(int2str(PTV1)) = 'Target1';
        organsMap(int2str(Brain)) = 'Brain';
        organsMap(int2str(BrainStem)) = 'BrainStem';
        organsMap(int2str(OpticNerves)) = 'OpticNerves';
        organsMap(int2str(Eyes)) = 'Eyes';
        organsMap(int2str(Parotids)) = 'Parotids';
        organsMap(int2str(Cord)) = 'Cord';
        organsMap(int2str(PTV2)) = 'Target2';
        organsMap(int2str(N2A)) = 'N2A';
        oMap = organsMap;
%     elseif patID == 5
%         caseParams.Skin = Skin;
%         caseParams.PTV1 = PTV1;
%         caseParams.Brain = Brain;
%         caseParams.BrainStem = BrainStem;
%         caseParams.OpticNerves = OpticNerves;
%         caseParams.Eyes = Eyes;
%         caseParams.Parotids = Parotids;
%         caseParams.Cord = Cord;
%         caseParams.PTV2 = PTV2;
%         caseParams.N2A_jxl = N2A_jxl;
%         caseParams.PTV3 = PTV3;
%         
%         organsMap(int2str(Skin)) = 'skin';
%         organsMap(int2str(PTV1)) = 'CTV 7060';
%         organsMap(int2str(Brain)) = 'Brain';
%         organsMap(int2str(BrainStem)) = 'BrainStem';
%         organsMap(int2str(OpticNerves)) = 'OpticNerves';
%         organsMap(int2str(Eyes)) = 'Eyes';
%         organsMap(int2str(Parotids)) = 'Parotids';
%         organsMap(int2str(Cord)) = 'Cord';
%         organsMap(int2str(PTV2)) = 'CTV 5040';
%         organsMap(int2str(N2A_jxl)) = 'N2A_jxl';
%         oMap = organsMap;
        
    elseif patID == 6
        caseParams.Skin = Skin;
        caseParams.PTV1 = PTV1;
        caseParams.Brain = Brain;
        caseParams.BrainStem = BrainStem;
        caseParams.OpticNerves = OpticNerves;
        caseParams.Eyes = Eyes;
        caseParams.LtParotids = LtParotids;
        caseParams.Cord = Cord;
        caseParams.PTV2 = PTV2;
        caseParams.RtParotids = RtParotids;
        caseParams.OpticChiasm = OpticChiasm;
        caseParams.Mandible = Mandible;
        caseParams.PTV3 = PTV3;
        
        
        organsMap(int2str(Skin)) = 'skin';
        organsMap(int2str(PTV3)) = 'Target3';
        organsMap(int2str(PTV1)) = 'Target1';
        organsMap(int2str(Brain)) = 'Brain';
        organsMap(int2str(BrainStem)) = 'BrainStem';
        organsMap(int2str(OpticNerves)) = 'OpticNerves';
        organsMap(int2str(Eyes)) = 'Eyes';
        organsMap(int2str(LtParotids)) = 'LtParotids';
        organsMap(int2str(Cord)) = 'Cord';
        organsMap(int2str(PTV2)) = 'Target2';
        organsMap(int2str(RtParotids)) = 'RtParotids';
        organsMap(int2str(Mandible)) = 'Mandible';
        organsMap(int2str(OpticChiasm)) = 'OpticChiasm';
        oMap = organsMap;
    elseif patID == 46
        caseParams.Brain = Brain;
        caseParams.BrainStem = BrainStem;
        caseParams.ExtAvoidance = ExtAvoidance;
        caseParams.Ext = Ext;
        caseParams.LtEye = LtEye;
        caseParams.Mandible = Mandible;
        caseParams.OpticChiasm = OpticChiasm;
        caseParams.OpticNerves = OpticNerves;
        caseParams.Parotids = Parotids;
        caseParams.RtEye = RtEye;
        caseParams.Skin = Skin;
        caseParams.Cord = Cord;
        caseParams.PTV1 = PTV1;
        caseParams.PTV2 = PTV2;
        caseParams.N3Jx1 = N3Jx1;
        caseParams.PTV3 = PTV3;
        
        
        organsMap(int2str(Brain)) = 'Brain';
        organsMap(int2str(BrainStem)) = 'BrainStem';
        organsMap(int2str(ExtAvoidance)) = 'ExtAvoidance';
        organsMap(int2str(Ext)) = 'Ext';
        organsMap(int2str(LtEye)) = 'LtEye';
        organsMap(int2str(Mandible)) = 'Mandible';
        organsMap(int2str(OpticChiasm)) = 'OpticChiasm';
        organsMap(int2str(OpticNerves)) = 'OpticNerves';
        organsMap(int2str(Parotids)) = 'Parotids';
        organsMap(int2str(RtEye)) = 'RtEye';
        organsMap(int2str(Skin)) = 'Skin';
        organsMap(int2str(Cord)) = 'Cord';
        organsMap(int2str(PTV1)) = 'PTV1';
        organsMap(int2str(PTV2)) = 'PTV2';
        organsMap(int2str(PTV3)) = 'PTV3';
        organsMap(int2str(N3Jx1)) = 'N3Jx1';
        
        oMap = organsMap;
    elseif patID == 8
        caseParams.Brain = Brain;
        caseParams.BrainStem = BrainStem;
        caseParams.Ext = Ext;
        caseParams.Mandible = Mandible;
        caseParams.LtLens = LtLens;
        caseParams.OpticChiasm = OpticChiasm;
        caseParams.OpticNerves = OpticNerves;
        caseParams.Eyes = Eyes;
        caseParams.Parotids = Parotids;
        caseParams.RtLens = RtLens;
        caseParams.Skin = Skin;
        caseParams.Cord = Cord;
        caseParams.PTV1 = PTV1;
        caseParams.PTV2 = PTV2;
        caseParams.PTV3 = PTV3;
        
        
        organsMap(int2str(Brain)) = 'Brain';
        organsMap(int2str(BrainStem)) = 'BrainStem';
        organsMap(int2str(Ext)) = 'Ext';
        organsMap(int2str(Mandible)) = 'Mandible';
        organsMap(int2str(LtLens)) = 'Left Lens';
        organsMap(int2str(OpticChiasm)) = 'OpticChiasm';
        organsMap(int2str(OpticNerves)) = 'OpticNerves';
        organsMap(int2str(Eyes)) = 'Eyes';
        organsMap(int2str(Parotids)) = 'Parotids';
        organsMap(int2str(RtLens)) = 'Right Lens';
        organsMap(int2str(Skin)) = 'Skin';
        organsMap(int2str(Cord)) = 'Cord';
        organsMap(int2str(PTV1)) = 'Target1';
        organsMap(int2str(PTV2)) = 'Target2';
        organsMap(int2str(PTV3)) = 'Target3';
        
        oMap = organsMap;
        
    elseif patID == 25
        caseParams.Brain = Brain;
        caseParams.BrainStem = BrainStem;
        caseParams.Ext = Ext;
        caseParams.LtEye = LtEye;
        caseParams.OpticNerves = OpticNerves;
        caseParams.Parotids = Parotids;
        caseParams.RtEye = RtEye;
        caseParams.Cord = Cord;
        caseParams.Skin = Skin;
        caseParams.PTV1 = PTV1;
        caseParams.PTV2 = PTV2;
        caseParams.PTV3 = PTV3;
        
        
        organsMap(int2str(Brain)) = 'Brain';
        organsMap(int2str(BrainStem)) = 'BrainStem';
        organsMap(int2str(Ext)) = 'Ext';
        organsMap(int2str(LtEye)) = 'LtEye';
        organsMap(int2str(OpticNerves)) = 'OpticNerves';
        organsMap(int2str(Parotids)) = 'Parotids';
        organsMap(int2str(RtEye)) = 'RtEye';
        organsMap(int2str(Cord)) = 'Cord';
        organsMap(int2str(Skin)) = 'Skin';
        organsMap(int2str(PTV1)) = 'PTV1';
        organsMap(int2str(PTV2)) = 'PTV2';
        organsMap(int2str(PTV3)) = 'PTV3';
        
        oMap = organsMap;
    elseif patID == 104
        caseParams.Brain = Brain;
        caseParams.BrainStem = BrainStem;
        caseParams.Ext = Ext;
        caseParams.LtParotid = LtParotid;
        caseParams.Mandible = Mandible;
        caseParams.OpticChiasm = OpticChiasm;
        caseParams.OpticNerves = OpticNerves;
        caseParams.Oralcavity = Oralcavity;
        caseParams.Eyes = Eyes;
        caseParams.RtParotid = RtParotid;
        caseParams.Skin = Skin;
        caseParams.Cord = Cord;
        caseParams.PTV1 = PTV1;
        caseParams.PTV2 = PTV2;
        caseParams.PTV3 = PTV3;
        
        
        organsMap(int2str(Brain)) = 'Brain';
        organsMap(int2str(BrainStem)) = 'BrainStem';
        organsMap(int2str(Ext)) = 'Ext';
        organsMap(int2str(LtParotid)) = 'LtParotid';
        organsMap(int2str(Mandible)) = 'Mandible';
        organsMap(int2str(OpticChiasm)) = 'OpticChiasm';
        organsMap(int2str(OpticNerves)) = 'OpticNerves';
        organsMap(int2str(Oralcavity)) = 'Oralcavity';
        organsMap(int2str(Eyes)) = 'Eyes';
        organsMap(int2str(RtParotid)) = 'RtParotid';
        organsMap(int2str(Skin)) = 'Skin';
        organsMap(int2str(Cord)) = 'Cord';
        organsMap(int2str(PTV1)) = 'CTV 7000';
        organsMap(int2str(PTV2)) = 'CTV 5400';
        organsMap(int2str(PTV3)) = 'PTV3';
        
        oMap = organsMap;
    elseif patID == 108
        caseParams.Brain = Brain;
        caseParams.BrainStem = BrainStem;
        caseParams.Ext = Ext;
        caseParams.LtParotid = LtParotid;
        caseParams.Mandible = Mandible;
        caseParams.OpticChiasm = OpticChiasm;
        caseParams.OpticNerves = OpticNerves;
        caseParams.Oralcavity = Oralcavity;
        caseParams.Eyes = Eyes;
        caseParams.RtParotid = RtParotid;
        caseParams.Skin = Skin;
        caseParams.Cord = Cord;
        caseParams.PTV1 = PTV1;
        caseParams.PTV2 = PTV2;
        caseParams.PTV3 = PTV3;
        
        
        organsMap(int2str(Brain)) = 'Brain';
        organsMap(int2str(BrainStem)) = 'BrainStem';
        organsMap(int2str(Ext)) = 'Ext';
        organsMap(int2str(LtParotid)) = 'LtParotid';
        organsMap(int2str(Mandible)) = 'Mandible';
        organsMap(int2str(OpticChiasm)) = 'OpticChiasm';
        organsMap(int2str(OpticNerves)) = 'OpticNerves';
        organsMap(int2str(Oralcavity)) = 'Oralcavity';
        organsMap(int2str(Eyes)) = 'Eyes';
        organsMap(int2str(RtParotid)) = 'RtParotid';
        organsMap(int2str(Skin)) = 'Skin';
        organsMap(int2str(Cord)) = 'Cord';
        organsMap(int2str(PTV1)) = 'CTV 7000';
        organsMap(int2str(PTV2)) = 'CTV 6000';
        organsMap(int2str(PTV3)) = 'PTV3';
        
        oMap = organsMap;
    elseif patID == 113
        caseParams.Brain = Brain;
        caseParams.BrainStem = BrainStem;
        caseParams.Ext = Ext;
        caseParams.LtParotid = LtParotid;
        caseParams.Mandible = Mandible;
        caseParams.OpticChiasm = OpticChiasm;
        caseParams.OpticNerves = OpticNerves;
        caseParams.Oralcavity = Oralcavity;
        caseParams.Eyes = Eyes;
        caseParams.RtParotid = RtParotid;
        caseParams.Skin = Skin;
        caseParams.Cord = Cord;
        caseParams.PTV1 = PTV1;
        caseParams.PTV2 = PTV2;
        caseParams.PTV3 = PTV3;
        
        
        organsMap(int2str(Brain)) = 'Brain';
        organsMap(int2str(BrainStem)) = 'BrainStem';
        organsMap(int2str(Ext)) = 'Ext';
        organsMap(int2str(LtParotid)) = 'LtParotid';
        organsMap(int2str(Mandible)) = 'Mandible';
        organsMap(int2str(OpticChiasm)) = 'OpticChiasm';
        organsMap(int2str(OpticNerves)) = 'OpticNerves';
        organsMap(int2str(Oralcavity)) = 'Oralcavity';
        organsMap(int2str(Eyes)) = 'Eyes';
        organsMap(int2str(RtParotid)) = 'RtParotid';
        organsMap(int2str(Skin)) = 'Skin';
        organsMap(int2str(Cord)) = 'Cord';
        organsMap(int2str(PTV1)) = 'CTV 7000';
        organsMap(int2str(PTV2)) = 'CTV 6000';
        organsMap(int2str(PTV3)) = 'PTV3';
        
        oMap = organsMap;
    elseif patID == 123
        caseParams.Brain = Brain;
        caseParams.BrainStem = BrainStem;
        caseParams.Ext = Ext;
        caseParams.LtParotid = LtParotid;
        caseParams.Mandible = Mandible;
        caseParams.OpticChiasm = OpticChiasm;
        caseParams.OpticNerves = OpticNerves;
        caseParams.Oralcavity = Oralcavity;
        caseParams.Eyes = Eyes;
        caseParams.RtParotid = RtParotid;
        caseParams.Skin = Skin;
        caseParams.Cord = Cord;
        caseParams.PTV1 = PTV1;
        caseParams.PTV2 = PTV2;
        caseParams.PTV3 = PTV3;
        
        
        organsMap(int2str(Brain)) = 'Brain';
        organsMap(int2str(BrainStem)) = 'BrainStem';
        organsMap(int2str(Ext)) = 'Ext';
        organsMap(int2str(LtParotid)) = 'LtParotid';
        organsMap(int2str(Mandible)) = 'Mandible';
        organsMap(int2str(OpticChiasm)) = 'OpticChiasm';
        organsMap(int2str(OpticNerves)) = 'OpticNerves';
        organsMap(int2str(Oralcavity)) = 'Oralcavity';
        organsMap(int2str(Eyes)) = 'Eyes';
        organsMap(int2str(RtParotid)) = 'RtParotid';
        organsMap(int2str(Skin)) = 'Skin';
        organsMap(int2str(Cord)) = 'Cord';
        organsMap(int2str(PTV1)) = 'CTV 7000';
        organsMap(int2str(PTV2)) = 'CTV 5600';
        organsMap(int2str(PTV3)) = 'PTV3';
        
        oMap = organsMap;
    elseif patID == 208
        caseParams.Brain = Brain;
        caseParams.BrainStem = BrainStem;
        caseParams.Ext = Ext;
        caseParams.ExtAvoidance = ExtAvoidance;
        caseParams.Larynx = Larynx;
        caseParams.LtParotid = LtParotid;
        caseParams.Mandible = Mandible;
        caseParams.MandiblePV = MandiblePV;
        caseParams.OpticChiasm = OpticChiasm;
        caseParams.OpticNerves = OpticNerves;
        caseParams.OpticNervesPV = OpticNervesPV;
        caseParams.OralCavity = OralCavity;
        caseParams.Eyes = Eyes;
        caseParams.RtParotid = RtParotid;
        caseParams.RtParotd = RtParotd;
        caseParams.Skin = Skin;
        caseParams.SkinPV = SkinPV;
        caseParams.Cord = Cord;
        caseParams.CordPV = CordPV;
        caseParams.PTV1 = PTV1;
        caseParams.PTV2 = PTV2;
        caseParams.PTV3 = PTV3;
        
        
        organsMap(int2str(Brain)) = 'Brain';
        organsMap(int2str(BrainStem)) = 'BrainStem';
        organsMap(int2str(Ext)) = 'Ext';
        organsMap(int2str(ExtAvoidance)) = 'ExtAvoidance';
        organsMap(int2str(Larynx)) = 'Larynx';
        organsMap(int2str(LtParotid)) = 'LtParotid';
        organsMap(int2str(Mandible)) = 'Mandible';
        organsMap(int2str(OpticChiasm)) = 'OpticChiasm';
        organsMap(int2str(OpticNerves)) = 'OpticNerves';
        organsMap(int2str(OralCavity)) = 'Oralcavity';
        organsMap(int2str(Eyes)) = 'Eyes';
        organsMap(int2str(RtParotid)) = 'RtParotid';
        organsMap(int2str(RtParotd)) = 'RtParotd';
        organsMap(int2str(Skin)) = 'Skin';
        organsMap(int2str(SkinPV)) = 'SkinPV';
        organsMap(int2str(Cord)) = 'Cord';
        organsMap(int2str(CordPV)) = 'CordPV';
        organsMap(int2str(PTV1)) = 'CTV 7000';
        organsMap(int2str(PTV2)) = 'CTV 6000';
        organsMap(int2str(PTV3)) = 'PTV3';
        
        oMap = organsMap;
    elseif patID == 30
        caseParams.Brain = Brain;
        caseParams.BrainStem = BrainStem;
        caseParams.LtEye = LtEye;
        caseParams.Mandible = Mandible;
        caseParams.OpticChiasm = OpticChiasm;
        caseParams.OpticNerves = OpticNerves;
        caseParams.Parotids = Parotids;
        caseParams.RtEye = RtEye;
        caseParams.Skin = Skin;
        caseParams.Cord = Cord;
        caseParams.PTV1 = PTV1;
        caseParams.PTV2 = PTV2;
        caseParams.PTV3 = PTV3;
        
        
        organsMap(int2str(Brain)) = 'Brain';
        organsMap(int2str(BrainStem)) = 'BrainStem';
        organsMap(int2str(LtEye)) = 'LtEye';
        organsMap(int2str(Mandible)) = 'Mandible';
        organsMap(int2str(OpticChiasm)) = 'OpticChiasm';
        organsMap(int2str(OpticNerves)) = 'OpticNerves';
        organsMap(int2str(Parotids)) = 'Parotids';
        organsMap(int2str(RtEye)) = 'RtEye';
        organsMap(int2str(Skin)) = 'Skin';
        organsMap(int2str(Cord)) = 'Cord';
        organsMap(int2str(PTV1)) = 'CTV 7030';
        organsMap(int2str(PTV2)) = 'CTV 6290';
        organsMap(int2str(PTV3)) = 'PTV3';
        
        oMap = organsMap;
    elseif patID == 79
        caseParams.Brain = Brain;
        caseParams.BrainStem = BrainStem;
        caseParams.LtParotid = LtParotid;
        caseParams.Mandible = Mandible;
        caseParams.OpticChiasm = OpticChiasm;
        caseParams.OpticNerves = OpticNerves;
        caseParams.Eyes = Eyes;
        caseParams.RtParotid = RtParotid;
        caseParams.Skin = Skin;
        caseParams.Cord = Cord;
        caseParams.PTV1 = PTV1;
        caseParams.PTV2 = PTV2;
        caseParams.PTV3 = PTV3;
        
        
        organsMap(int2str(Brain)) = 'Brain';
        organsMap(int2str(BrainStem)) = 'BrainStem';
        organsMap(int2str(LtParotid)) = 'LtParotid';
        organsMap(int2str(Mandible)) = 'Mandible';
        organsMap(int2str(OpticChiasm)) = 'OpticChiasm';
        organsMap(int2str(OpticNerves)) = 'OpticNerves';
        organsMap(int2str(Eyes)) = 'Eyes';
        organsMap(int2str(RtParotid)) = 'RtParotid';
        organsMap(int2str(Skin)) = 'Skin';
        organsMap(int2str(Cord)) = 'Cord';
        organsMap(int2str(PTV1)) = 'CTV 7030';
        organsMap(int2str(PTV2)) = 'CTV 6290';
        organsMap(int2str(PTV3)) = 'PTV3';
        
        oMap = organsMap;
    elseif patID == 17
        caseParams.BrainStem = BrainStem;
        caseParams.Ext = Ext;
        caseParams.LtEye = LtEye;
        caseParams.OpticNerves = OpticNerves;
        caseParams.Parotids = Parotids;
        caseParams.RtEye = RtEye;
        caseParams.Skin = Skin;
        caseParams.Cord = Cord;
        caseParams.PTV1 = PTV1;
        caseParams.PTV2 = PTV2;
        caseParams.PTV3 = PTV3;
              
        
        organsMap(int2str(BrainStem)) = 'BrainStem';
        organsMap(int2str(Ext)) = 'Ext';
        organsMap(int2str(LtEye)) = 'LtEye';
        organsMap(int2str(OpticNerves)) = 'OpticNerves';
        organsMap(int2str(Parotids)) = 'Parotids';
        organsMap(int2str(RtEye)) = 'RtEye';
        organsMap(int2str(Skin)) = 'Skin';
        organsMap(int2str(Cord)) = 'Cord';
        organsMap(int2str(PTV1)) = 'CTV 7030';
        organsMap(int2str(PTV2)) = 'CTV 5940';
        organsMap(int2str(PTV3)) = 'PTV3';
        
        oMap = organsMap;
    else
        prostrate_param
        oMap = organsMap;
    end
    caseParams.hotspot = hotspot;
    %caseParams.anchor = anchor;
    %caseParams.inHotspot = inHotspot;
    
    if exist('rectumSlices') %for type N
        caseParams.rectumSlices = rectumSlices;
    else
        caseParams.rectumSlices = [];
        warning('No rectum slices found in init_pat file.  Not needed for type ''orig''.');
    end
    
    if exist('prescriptionDose') %for type P
        caseParams.prescriptionDose = prescriptionDose;
    else
        %caseParams.prescriptionDose = [];
        warning('No prescription dose found in init_pat file.  Only needed for type ''P''.');
    end
catch
    error('One or more of the case-specific parameters necessary to run prioritized prescription optimization have not been defined.');
end
