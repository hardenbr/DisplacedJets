# sum quantities

python plot1D.py -f calolist.py -v "jetIPSigLogSum2D" -c "1" -t jets --xmin -20 --xmax 60 --nbins 100 --xlabel "Jet Sum over Log(IP_{sig})" -o test.root --genmatch 
python plot1D.py -f calolist.py -v "jetDistLogSum" -c "1" -t jets --xmin 0 --xmax 200 --nbins 50 --xlabel "Jet Sum over Log(Dist)" -o test.root --genmatch 
python plot1D.py -f calolist.py -v "jetIPSigLogSum3D" -c "1" -t jets --xmin -100 --xmax 200 --nbins 50 --xlabel "Jet Sum over Log(IP_{sig})" -o test.root --genmatch 

# median
python plot1D.py -f calolist.py -v "jetMedianIPSig2D" -c "1" -t jets --xmin 0 --xmax 1 --nbins 50  --xlabel "Jet Median 2D IP_{sig}" -o test.root --genmatch 
python plot1D.py -f calolist.py -v "jetMedianIPLogSig2D" -c "1" -t jets --xmin -10 --xmax 10 --nbins 40 --xlabel "Jet Median IP Log Sig 2D" -o test.root --genmatch 

#variance
python plot1D.py -f calolist.py -v "jetVarianceIPSig2D" -c "1" -t jets --xmin 0 --xmax 1000 --nbins 100 --xlabel "IP_{sig} Variance" -o test.root --genmatch 
python plot1D.py -f calolist.py -v "jetVarianceJetDist" -c "1" -t jets --xmin 0 --xmax 1000 --nbins 100 --xlabel "Track Jet Distance Variance" -o test.root --genmatch 
python plot1D.py -f calolist.py -v "jetVarianceIP2D" -c "1" -t jets --xmin 0 --xmax 1000 --nbins 100 --xlabel "IP2D Variance" -o test.root --genmatch
python plot1D.py -f calolist.py -v "jetVarianceIPLogSig2D" -c "1" -t jets --xmin -100 --xmax 100 --nbins 100 --xlabel "IP Log Sig 2D Variance" -o test.root --genmatch 

### mean 
python plot1D.py -f calolist.py -v "jetMeanIPSig2D" -c "1" -t jets --xmin 0 --xmax 1000 --nbins 100 --xlabel "IP_{sig} Mean" -o test.root --genmatch 
python plot1D.py -f calolist.py -v "jetMeanJetDist" -c "1" -t jets --xmin 0 --xmax 1000 --nbins 100 --xlabel "Jet Track Distance Mean" -o test.root --genmatch 
python plot1D.py -f calolist.py -v "jetMeanIPLogSig2D" -c "1" -t jets --xmin -10 --xmax 10 --nbins 40 --xlabel "Jet Mean IP Log Sig 2D" -o test.root --genmatch 

### secondary vertex
python plot1D.py -f calolist.py -v "jetSvMass" -c "1" -t jets --xmin 0 --xmax 20 --nbins 100 --xlabel "SV Mass " -o test.root --genmatch 
python plot1D.py -f calolist.py -v "jetSvNTrack" -c "1" -t jets --xmin -.5 --xmax 19.5 --nbins 20 --xlabel "SV Track Multiplicity" -o test.root --genmatch 
python plot1D.py -f calolist.py -v "jetSvLxySig" -c "1" -t jets --xmin 0 --xmax 200 --nbins 100 --xlabel "SV L_{xy} Significance" -o test.root --genmatch 
python plot1D.py -f calolist.py -v "jetSvLxy" -c "1" -t jets --xmin 0 --xmax 20 --nbins 100 --xlabel "SV L_{xy}" -o test.root --genmatch 
python plot1D.py -f calolist.py -v "jetNSv" -c "1" -t jets --xmin -.5 --xmax 9.5 --nbins 10 --xlabel "N Reco SV" -o test.root --genmatch 

### IVF MATCHED VERTEX
python plot1D.py -f calolist.py -v "jetIVFMass" -c "1" -t jets --xmin 0 --xmax 20 --nbins 100 --xlabel "IVF Mass " -o test.root 
python plot1D.py -f calolist.py -v "jetIVFNTrack" -c "1" -t jets --xmin -.5 --xmax 19.5 --nbins 20 --xlabel "IVF Track Multiplicity" -o test.root
python plot1D.py -f calolist.py -v "jetIVFLxySig" -c "1" -t jets --xmin 0 --xmax 200 --nbins 100 --xlabel "IVF L_{xy} Significance" -o test.root
python plot1D.py -f calolist.py -v "jetIVFLxy" -c "1" -t jets --xmin 0 --xmax 20 --nbins 100 --xlabel "IVF L_{xy}" -o test.root
python plot1D.py -f calolist.py -v "jetIVFMatchingScore" -c "1" -t jets --xmin 0 --xmax 20 --nbins 100 --xlabel "IVF Matching Score" -o test.root 

## calo jet information
python plot1D.py -f calolist.py -v "caloJetHfrac" -c "1" -t jets --xmin 0 --xmax 1 --nbins 100 --xlabel "Hadronic Energy Fraction" -o test.root --genmatch 
python plot1D.py -f calolist.py -v "caloJetPt" -c "1" -t jets --xmin 40 --xmax 700 --nbins 100 --xlabel "Calo Jet p_{t} [GeV]" -o test.root --genmatch 
python plot1D.py -f calolist.py -v "caloJetEta" -c "1" -t jets --xmin -3 --xmax 3 --nbins 100 --xlabel "Calo Jet #eta" -o test.root --genmatch 
python plot1D.py -f calolist.py -v "caloJetPhi" -c "1" -t jets --xmin -3 --xmax 3 --nbins 100 --xlabel "Calo Jet #phi" -o test.root --genmatch 
python plot1D.py -f calolist.py -v "caloJetTowerArea" -c "1" -t jets --xmin 0 --xmax 1 --nbins 100 --xlabel "Tower Area" -o test.root --genmatch 

#signal 2d
python plot2D.py -f sig2D.list -v "jetSvNTrack:jetMedianIPSig2D" -c "1" -t jets --xmin 0 --xmax 5 --nxbins 25  --xlabel "Jet Median 2D IP Sig." --ymin 0.5 --ymax 15.5 --nybins 15 --ylabel "N SV Tracks" -o test.root --genmatch
python plot2D.py -f sig2D.list -v "jetIPSigLogSum2D:jetMedianIPSig2D" -c "1" -t jets --xmin 0 --xmax 4 --nxbins 25  --xlabel "Jet Median 2D IP_{sig}" --ymin -10 --ymax 100 --nybins 25 --ylabel "Jet Sum over Log(IP_{sig})" -o test.root --genmatch
python plot2D.py -f sig2D.list -v "jetIPSigLogSum2D:jetSvLxySig" -c "1" -t jets --xmin 0 --xmax 300 --nxbins 30  --xlabel "Jet Secondary Vertex L_{xy} Sig." --ymin -10 --ymax 100 --nybins 25 --ylabel "Jet Sum over Log(IP_{sig})" -o test.root --genmatch
python plot2D.py -f sig2D.list -v "jetIPSigLogSum2D:jetSvNTrack" -c "1" -t jets --xmin 0 --xmax 15 --nxbins 15  --xlabel "SV Track Multiplicty" --ymin -10 --ymax 100 --nybins 25 --ylabel "Jet Sum over Log(IP_{sig})" -o test.root --genmatch
python plot2D.py -f sig2D.list -v "jetSvNTrack:jetSvLxySig" -c "1" -t jets --xmin 0 --xmax 300 --nxbins 30  --xlabel "Jet Secondary Vertex L_{xy} Sig." --ymin 0 --ymax 20 --nybins 20 --ylabel "SV Track Multiplicity" -o test.root --genmatch
python plot2D.py -f sig2D.list -v "jetSvNTrack:jetSvMass" -c "1" -t jets --xmin 0 --xmax 100 --nxbins 30  --xlabel "Jet Secondary Vertex Mass [GeV]" --ymin 0 --ymax 20 --nybins 20 --ylabel "SV Track Multiplicity" -o test.root --genmatch
python plot2D.py -f sig2D.list -v "jetIPLogSum2D:jetIPSigLogSum2D" -c "1" -t jets --xmin -20 --xmax 100 --nxbins 60  --xlabel "Jet Sum over  Log(2DIP_{sig}})" --ymin -300 --ymax 50 --nybins 60 --ylabel "Jet Sum over Log(2DIP)" -o test.root --genmatch

#background 2d
#IP SIG LOG SUM VS XXXX
python plot2D.py -f qcd2D.list -v "jetIPSigLogSum2D:jetMedianIPSig2D" -c "1" -t jets --xmin 0 --xmax 4 --nxbins 25  --xlabel "Jet Median 2D IP_{sig}" --ymin -10 --ymax 100 --nybins 25 --ylabel "Jet Sum over Log(IP_{sig})" -o test.root --genmatch
python plot2D.py -f qcd2D.list -v "jetIPSigLogSum2D:jetIVFLxySig" -c "1" -t jets --xmin 0 --xmax 500 --nxbins 30  --xlabel "Jet IVF L_{xy} Sig" --ymin -10 --ymax 100 --nybins 25 --ylabel "Jet Sum over Log(IP_{sig})" -o test.root --genmatch
python plot2D.py -f qcd2D.list -v "jetIPSigLogSum2D:jetIVFNTrack" -c "1" -t jets --xmin 0 --xmax 500 --nxbins 30  --xlabel "IVF N Tracks" --ymin -10 --ymax 100 --nybins 25 --ylabel "Jet Sum over Log(IP_{sig})" -o test.root --genmatch
python plot2D.py -f qcd2D.list -v "jetIPSigLogSum2D:jetDistLogSum" -c "1" -t jets --xmin 0 --xmax 200 --nxbins 30  --xlabel "Jet Dist Log Sum" --ymin -10 --ymax 100 --nybins 25 --ylabel "Jet Sum over Log(IP_{sig})" -o test.root --genmatch
python plot2D.py -f qcd2D.list -v "jetIPSigLogSum2D:jetMeanIPSig2D" -c "1" -t jets --xmin 0 --xmax 1000 --nxbins 30  --xlabel "Jet Mean IP Sig 2D" --ymin -10 --ymax 100 --nybins 25 --ylabel "Jet Sum over Log(IP_{sig})" -o test.root --genmatch
python plot2D.py -f qcd2D.list -v "jetIPSigLogSum2D:jetVarianceIPSig2D" -c "1" -t jets --xmin 0 --xmax 200000 --nxbins 30  --xlabel "Jet Variance IP Sig 2D" --ymin -50 --ymax 250 --nybins 25 --ylabel "Jet Sum Log(IP_{sig})" -o test.root --genmatch
python plot2D.py -f qcd2D.list -v "jetIPSigLogSum2D:jetMedianIPLogSig2D" -c "1" -t jets --xmin -5 --xmax 10 --nxbins 30  --xlabel "Jet Median IP Log Sig 2D" --ymin -50 --ymax 250 --nybins 25 --ylabel "Jet Sum Log(IP_{sig})" -o test.root --genmatch
#IP SIG MEDIAN VS XX
 python plot2D.py -f qcd2D.list -v "jetMedianIPSig2D:jetVarianceIPSig2D" -c "1" -t jets --xmin 0 --xmax 200000 --nxbins 30  --xlabel "Jet Variance IP Sig 2D" --ymin 0 --ymax 10 --nybins 25 --ylabel "Jet Median IP Sig 2D" -o test.root --genmatch
 python plot2D.py -f qcd2D.list -v "jetMedianIPSig2D:jetIVFLxySig" -c "1" -t jets --xmin 0 --xmax 3000 --nxbins 30  --xlabel "Jet IVF Lxy Sig" --ymin 0 --ymax 10 --nybins 25 --ylabel "Jet Median IP Sig 2D" -o test.root --ivfvtxmatch
 python plot2D.py -f qcd2D.list -v "jetMedianIPSig2D:jetMedianIPLogSig2D" -c "1" -t jets --xmin -5 --xmax 10 --nxbins 30  --xlabel "Jet Median IP Log Sig 2D" --ymin 0 --ymax 10 --nybins 25 --ylabel "Jet Median IP Sig 2D" -o test.root --ivfvtxmatch
python plot2D.py -f qcd2D.list -v "jetMedianIPLogSig2D:jetIVFLxySig" -c "1" -t jets --xmin 0 --xmax 3000 --nxbins 30  --xlabel "Jet IVF Lxy Sig" --ymin -5 --ymax 10 --nybins 25 --ylabel "Jet Median IP Log Sig 2D" -o test.root --ivfvtxmatch
python plot2D.py -f qcd2D.list -v "jetMedianIPLogSig2D:jetIVFLxySig" -c "1" -t jets --xmin 0 --xmax 3000 --nxbins 30  --xlabel "Jet SV Lxy Sig" --ymin -5 --ymax 10 --nybins 25 --ylabel "Jet Median IP Log Sig 2D" -o test.root --ivfvtxmatch
#IVF LXY SIG
python plot2D.py -f qcd2D.list -v "jetIVFNTrack:jetIVFLxySig" -c "1" -t jets --xmin 0 --xmax 3000 --nxbins 30  --xlabel "Jet SV Lxy Sig" --ymin 0 --ymax 60 --nybins 60 --ylabel "Jet IVF N Tracks" -o test.root
python plot2D.py -f qcd2D.list -v "jetIVFNTrack:log(jetIVFLxySig)" -c "1" -t jets --xmin -2 --xmax 10 --nxbins 30  --xlabel "Jet SV Log Lxy Sig" --ymin 0 --ymax 60 --nybins 60 --ylabel "Jet IVF N Tracks" -o test.root
#IVF N TRACKS
python plot2D.py -f qcd2D.list -v "jetIVFNTrack:jetMedianIPLogSig2D" -c "1" -t jets --xmin -2 --xmax 8 --nxbins 30  --xlabel "Jet Median IP Log Sig2D" --ymin 0 --ymax 70 --nybins 30 --ylabel "Jet IVF N Tracks" -o test.root

# MISC CHECKS
python plot2D.py -f qcd2D.list -v "jetVarianceIPSig2D:jetMeanIPSig2D" -c "1" -t jets --xmin 0 --xmax 1000 --nxbins 30  --xlabel "Jet Mean IP Sig 2D" --ymin 0 --ymax 1000 --nybins 25 --ylabel "Jet Variance Log(IP_{sig})" -o test.root --genmatch

### PLOT TAGGING VARIABLES
python plot1D.py -f tags.list -v "simple" -c "simpleID==1" -t jets --xmin -.001 --xmax .0015 --nbins 100 --xlabel "Linear Fischer Discriminator" -o test.root --genmatch
python tag_analyzer.py -s signal.root -b qcd.root -t jets -o tags.root --gidout GID.txt 

### PLOT GENERATOR TREE
#generator level vertices restricted to the first event 
python plot2D.py -f vtx2D.list -v "genPartVY:genPartVX" -c "evNum==1" -t gen --xmin -5 --xmax 5 --nxbins 50  --xlabel "Gen Particle Vertex X [cm]" --ymin -5 --ymax 5 --nybins 50 --ylabel "Gen Particle Vertex Y [cm]" -o test.root
# simulated vertices restricted to vertices from the primary generator process
python plot2D.py -f vtx2D.list -v "simVtxY:simVtxX" -c "evNum==1 && simVtxProcType==0" -t gen --xmin -5 --xmax 5 --nxbins 50  --xlabel "Sim Generator Vertex X [cm]" --ymin -5 --ymax 5 --nybins 50 --ylabel "Sim Gnenerator Vertex Y [cm]" -o test.root

### PLOT VERTEX TREE
python plot2D.py -f vtx2D.list -v "vtxIncSecY:vtxIncSecX" -c "evNum==1" -t vtx --xmin -5 --xmax 5 --nxbins 50  --xlabel "Inclusive Secondary Vertex X [cm]" --ymin -5 --ymax 5 --nybins 50 --ylabel "Inclusive Secondary Vertex Y [cm]" -o test.root
python plot2D.py -f vtx2D.list -v "vtxIncSecY:vtxIncSecX" -c "vtxIncSecGenMatchMetric * (vtxIncSecGenMatched > 0 && evNum==1) " -t vtx --xmin -5 --xmax 5 --nxbins 50  --xlabel "Inclusive Secondary Vertex X [cm]" --ymin -5 --ymax 5 --nybins 50 --ylabel "Inclusive Secondary Vertex Y [cm]" -o test.root
python plot2D.py -f vtx2D.list -v "jetSvX:jetSvY" -c "evNum==1" -t jets --xmin -5 --xmax 5 --nxbins 50  --xlabel "AVF Vertex X [cm]" --ymin -5 --ymax 5 --nybins 50 --ylabel "AVF Vertex Y [cm]" -o test.root
