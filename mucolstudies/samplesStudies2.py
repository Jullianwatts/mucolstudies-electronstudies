import math
import glob
import ROOT
import os
from podio import root_io
exec(open("./plotHelper.py").read())
ROOT.gROOT.SetBatch()
PLOT_DIR = "/scratch/jwatts/mucol/mucolstudies/plots2026/endcap"
os.makedirs(PLOT_DIR, exist_ok=True)
samples = sorted(glob.glob("/scratch/jwatts/mucol/v3/reco/electronGun_0_50/*_reco_*.edm4hep.root"))
#samples = sorted(glob.glob("/scratch/jwatts/mucol/v3/reco/electronGun_50_250/*_reco_*.edm4hep.root"))
#samples = sorted(glob.glob("/scratch/jwatts/mucol/v3/reco/electronGun_250_1000/*_reco_*.edm4hep.root"))
#samples = sorted(glob.glob("/scratch/jwatts/mucol/v3/reco/electronGun_1000_5000/*_reco_*.edm4hep.root"))
files = {"electronGun_0_50": samples}

B_FIELD = 5
M_ELECTRON = 0.000511
DR_CONE = 0.2
ETA_SLICES = [(0.0, 0.5), (0.5, 1.0), (1.0, 1.5), (1.5, 2.0), (2.0, 2.4)]
FWD_LO = 1.5            # the region we want to understand
FWD_HI = 2.4
CEN_HI = 1.0            # central region to compare it against
N_PRINT_FWD = 25        # detailed dump of this many forward events

coll_candidates = {
    "mcp":      ["MCParticles", "MCParticle"],
    "trk":      ["SiTracks", "SelectedTracks", "SiTracks_Refitted"],
    "cluster":  ["PandoraClusters"],
    "pfo":      ["PandoraPFOs"],
}
coll_names = {}

def resolveCollections(event):
    available = list(event.getAvailableCollections())
    for key in coll_candidates:
        coll_names[key] = None
        for cand in coll_candidates[key]:
            if cand in available:
                coll_names[key] = cand
                break
        print(f"Using {key:<8} -> {coll_names[key]}")
    print()
    return

def getCollection(event, key):
    name = coll_names.get(key)
    if name is None: return []
    return event.get(name)

def getVec(v):
    return (v.x, v.y, v.z)

def getMCTLV(mcp):
    px, py, pz = getVec(mcp.getMomentum())
    m = mcp.getMass()
    E = math.sqrt(px**2 + py**2 + pz**2 + m**2)
    tlv = ROOT.TLorentzVector()
    tlv.SetPxPyPzE(px, py, pz, E)
    return tlv

def getRecoTLV(pfo):
    px, py, pz = getVec(pfo.getMomentum())
    E = pfo.getEnergy()
    tlv = ROOT.TLorentzVector()
    tlv.SetPxPyPzE(px, py, pz, E)
    return tlv

def getPDG(pfo):
    if hasattr(pfo, "getPDG"): return pfo.getPDG()
    return pfo.getType()

class TrackParams:
    # edm4hep keeps the helix parameters in the TrackState instead of on the
    # track itself, so wrap them to look like an LCIO track and hand it to
    # getPt/getP/getTrackTLV from plotHelper.py
    def __init__(self, trk):
        self.ts = None
        states = trk.getTrackStates()
        if len(states) > 0:
            self.ts = states[0]
            for ts in states:
                if ts.location == 1:   # 1 = AtIP
                    self.ts = ts
                    break
    def getOmega(self):     return self.ts.omega
    def getTanLambda(self): return self.ts.tanLambda
    def getPhi(self):       return self.ts.phi
    def getD0(self):        return self.ts.D0
    def getZ0(self):        return self.ts.Z0

def getTrackTLV_edm4hep(trk, m = M_ELECTRON, b_field = B_FIELD):
    params = TrackParams(trk)
    if params.ts is None or params.ts.omega == 0: return ROOT.TLorentzVector()
    return getTrackTLV(params, m, b_field)

def getClusterTLV_edm4hep(cluster):
    x, y, z = getVec(cluster.getPosition())
    E = cluster.getEnergy()
    r = math.sqrt(x**2 + y**2 + z**2)
    if r == 0: return ROOT.TLorentzVector()
    tlv = ROOT.TLorentzVector()
    tlv.SetPxPyPzE(E*x/r, E*y/r, E*z/r, E)
    return tlv

slice_labels = [f"{lo:.1f} - {hi:.1f}" for lo, hi in ETA_SLICES]
slice_keys = ["n", "trk", "clu", "pfo", "pfo_el", "pfo_ch", "trkclu", "sum_nclu", "sum_Eratio", "n_Eratio", "sum_dr_clu"]
slices = {}
for lab in slice_labels:
    slices[lab] = {k: 0 for k in slice_keys}

def getSliceLabel(eta):
    for i, (lo, hi) in enumerate(ETA_SLICES):
        if lo <= abs(eta) < hi: return slice_labels[i]
    return None

h = {}
# every reco object in the file, finely binned so the acceptance edge is visible
h["trk_eta_all"] = ROOT.TH1F("trk_eta_all", "", 120, -3, 3)
h["clu_eta_all"] = ROOT.TH1F("clu_eta_all", "", 120, -3, 3)
h["pfo_eta_all"] = ROOT.TH1F("pfo_eta_all", "", 120, -3, 3)
h["mcp_eta_all"] = ROOT.TH1F("mcp_eta_all", "", 120, -3, 3)
# truth electrons in the forward region, split by what was found
h["fwd_mcp_eta"]    = ROOT.TH1F("fwd_mcp_eta",    "", 18, 1.5, 2.4)
h["fwd_mcp_eta_trk"] = ROOT.TH1F("fwd_mcp_eta_trk", "", 18, 1.5, 2.4)
h["fwd_mcp_eta_clu"] = ROOT.TH1F("fwd_mcp_eta_clu", "", 18, 1.5, 2.4)
h["fwd_mcp_eta_pfo"] = ROOT.TH1F("fwd_mcp_eta_pfo", "", 18, 1.5, 2.4)
h["fwd_mcp_pt"]     = ROOT.TH1F("fwd_mcp_pt",     "", 50, 0, 50)
h["fwd_mcp_pt_trk"] = ROOT.TH1F("fwd_mcp_pt_trk", "", 50, 0, 50)
# central vs forward comparisons
h["cen_dr_trk"] = ROOT.TH1F("cen_dr_trk", "", 100, 0, 2)
h["fwd_dr_trk"] = ROOT.TH1F("fwd_dr_trk", "", 100, 0, 2)
h["cen_dr_clu"] = ROOT.TH1F("cen_dr_clu", "", 100, 0, 2)
h["fwd_dr_clu"] = ROOT.TH1F("fwd_dr_clu", "", 100, 0, 2)
h["cen_dr_pfo"] = ROOT.TH1F("cen_dr_pfo", "", 100, 0, 2)
h["fwd_dr_pfo"] = ROOT.TH1F("fwd_dr_pfo", "", 100, 0, 2)
h["cen_Eratio"] = ROOT.TH1F("cen_Eratio", "", 100, 0, 2)
h["fwd_Eratio"] = ROOT.TH1F("fwd_Eratio", "", 100, 0, 2)
h["cen_nclu"]   = ROOT.TH1F("cen_nclu",   "", 11, -0.5, 10.5)
h["fwd_nclu"]   = ROOT.TH1F("fwd_nclu",   "", 11, -0.5, 10.5)
h["cen_clu_E"]  = ROOT.TH1F("cen_clu_E",  "", 60, 0, 60)
h["fwd_clu_E"]  = ROOT.TH1F("fwd_clu_E",  "", 60, 0, 60)

pdg_counts = {"central": {}, "forward": {}}
n_events = 0
n_fwd_printed = 0
resolved = False

print(f"\n--- First {N_PRINT_FWD} events with a truth electron in {FWD_LO} < |#eta| < {FWD_HI} ---")
print(f"{'Event':<6} | {'eta':<6} | {'pT':<7} | {'nTrk':<4} | {'drTrk':<7} | {'nClu':<4} | "
      f"{'cluE':<7} | {'cluEta':<7} | {'drClu':<7} | {'PFO PDGs'}")

for s in files:
    for f in files[s]:
        reader = root_io.Reader(f)
        for event in reader.get("events"):
            if not resolved:
                print()
                resolveCollections(event)
                resolved = True
            n_events += 1
            mcps = getCollection(event, "mcp")
            trks = getCollection(event, "trk")
            clusters = getCollection(event, "cluster")
            pfos = getCollection(event, "pfo")

            mcp_electrons = []
            for mcp in mcps:
                if mcp.getGeneratorStatus() != 1 or abs(mcp.getPDG()) != 11: continue
                tlv = getMCTLV(mcp)
                mcp_electrons.append(tlv)
                h["mcp_eta_all"].Fill(tlv.Eta())

            trk_tlvs = []
            for t in trks:
                tlv = getTrackTLV_edm4hep(t)
                trk_tlvs.append(tlv)
                h["trk_eta_all"].Fill(tlv.Eta())

            clu_tlvs = []
            for c in clusters:
                tlv = getClusterTLV_edm4hep(c)
                clu_tlvs.append(tlv)
                h["clu_eta_all"].Fill(tlv.Eta())

            pfo_tlvs_all = []
            pfo_tlvs_el = []
            pfo_tlvs_ch = []
            pfo_pdgs = []
            for p in pfos:
                tlv = getRecoTLV(p)
                pdg = abs(getPDG(p))
                pfo_tlvs_all.append(tlv)
                pfo_pdgs.append(pdg)
                if pdg == 11: pfo_tlvs_el.append(tlv)
                if pdg in (11, 13, 211): pfo_tlvs_ch.append(tlv)
                h["pfo_eta_all"].Fill(tlv.Eta())

            for mcp_el in mcp_electrons:
                m_eta = mcp_el.Eta()
                m_pt = mcp_el.Perp()
                lab = getSliceLabel(m_eta)

                drs_trk = sorted([(mcp_el.DeltaR(t), t) for t in trk_tlvs], key=lambda x: x[0])
                drs_clu = sorted([(mcp_el.DeltaR(c), c) for c in clu_tlvs], key=lambda x: x[0])
                drs_pfo = sorted([(mcp_el.DeltaR(p), p) for p in pfo_tlvs_all], key=lambda x: x[0])
                drs_el  = sorted([(mcp_el.DeltaR(p), p) for p in pfo_tlvs_el], key=lambda x: x[0])
                drs_ch  = sorted([(mcp_el.DeltaR(p), p) for p in pfo_tlvs_ch], key=lambda x: x[0])

                has_trk = len(drs_trk) > 0 and drs_trk[0][0] < DR_CONE
                has_clu = len(drs_clu) > 0 and drs_clu[0][0] < DR_CONE
                has_pfo = len(drs_pfo) > 0 and drs_pfo[0][0] < DR_CONE
                has_el  = len(drs_el)  > 0 and drs_el[0][0]  < DR_CONE
                has_ch  = len(drs_ch)  > 0 and drs_ch[0][0]  < DR_CONE
                n_clu_cone = len([1 for dr, c in drs_clu if dr < DR_CONE])
                E_cone = sum([c.E() for dr, c in drs_clu if dr < DR_CONE])

                if lab is not None:
                    slices[lab]["n"] += 1
                    if has_trk: slices[lab]["trk"] += 1
                    if has_clu: slices[lab]["clu"] += 1
                    if has_pfo: slices[lab]["pfo"] += 1
                    if has_el:  slices[lab]["pfo_el"] += 1
                    if has_ch:  slices[lab]["pfo_ch"] += 1
                    if has_trk and has_clu: slices[lab]["trkclu"] += 1
                    slices[lab]["sum_nclu"] += n_clu_cone
                    if len(drs_clu) > 0: slices[lab]["sum_dr_clu"] += drs_clu[0][0]
                    if has_clu and mcp_el.E() > 0:
                        slices[lab]["sum_Eratio"] += E_cone/mcp_el.E()
                        slices[lab]["n_Eratio"] += 1

                is_fwd = FWD_LO <= abs(m_eta) < FWD_HI
                is_cen = abs(m_eta) < CEN_HI
                region = "forward" if is_fwd else ("central" if is_cen else None)
                if region is not None:
                    for pdg in pfo_pdgs:
                        pdg_counts[region][pdg] = pdg_counts[region].get(pdg, 0) + 1
                    pre = "fwd" if is_fwd else "cen"
                    if len(drs_trk) > 0: h[pre+"_dr_trk"].Fill(drs_trk[0][0])
                    if len(drs_clu) > 0: h[pre+"_dr_clu"].Fill(drs_clu[0][0])
                    if len(drs_pfo) > 0: h[pre+"_dr_pfo"].Fill(drs_pfo[0][0])
                    h[pre+"_nclu"].Fill(n_clu_cone)
                    if has_clu:
                        h[pre+"_Eratio"].Fill(E_cone/mcp_el.E())
                        h[pre+"_clu_E"].Fill(drs_clu[0][1].E())

                if is_fwd:
                    h["fwd_mcp_eta"].Fill(abs(m_eta))
                    h["fwd_mcp_pt"].Fill(m_pt)
                    if has_trk:
                        h["fwd_mcp_eta_trk"].Fill(abs(m_eta))
                        h["fwd_mcp_pt_trk"].Fill(m_pt)
                    if has_clu: h["fwd_mcp_eta_clu"].Fill(abs(m_eta))
                    if has_pfo: h["fwd_mcp_eta_pfo"].Fill(abs(m_eta))

                    if n_fwd_printed < N_PRINT_FWD:
                        n_fwd_printed += 1
                        dr_trk_s = f"{drs_trk[0][0]:.3f}" if len(drs_trk) > 0 else "-"
                        dr_clu_s = f"{drs_clu[0][0]:.3f}" if len(drs_clu) > 0 else "-"
                        clu_E_s  = f"{drs_clu[0][1].E():.2f}" if len(drs_clu) > 0 else "-"
                        clu_eta_s = f"{drs_clu[0][1].Eta():.2f}" if len(drs_clu) > 0 else "-"
                        print(f"{n_events:<6} | {m_eta:<6.2f} | {m_pt:<7.2f} | {len(trk_tlvs):<4} | {dr_trk_s:<7} | "
                              f"{len(clu_tlvs):<4} | {clu_E_s:<7} | {clu_eta_s:<7} | {dr_clu_s:<7} | {sorted(pfo_pdgs)}")

atltext = ["Muon Collider", "Simulation, no BIB"]

plotHistograms({"truth electron": h["mcp_eta_all"], "track": h["trk_eta_all"],
                "cluster": h["clu_eta_all"], "PFO": h["pfo_eta_all"]},
               os.path.join(PLOT_DIR, "eta_acceptance.png"), "#eta", "Objects", atltext=atltext)
plotHistograms({"all truth": h["fwd_mcp_eta"], "with track": h["fwd_mcp_eta_trk"],
                "with cluster": h["fwd_mcp_eta_clu"], "with PFO": h["fwd_mcp_eta_pfo"]},
               os.path.join(PLOT_DIR, "forward_eta_matched.png"), "|#eta|", "Truth electrons", atltext=atltext)
plotHistograms({"all truth": h["fwd_mcp_pt"], "with track": h["fwd_mcp_pt_trk"]},
               os.path.join(PLOT_DIR, "forward_pt_matched.png"), "p_{T} [GeV]", "Truth electrons", atltext=atltext)
plotHistograms({f"|#eta| < {CEN_HI}": h["cen_dr_trk"], f"{FWD_LO} < |#eta| < {FWD_HI}": h["fwd_dr_trk"]},
               os.path.join(PLOT_DIR, "dr_track_by_region.png"), "#DeltaR(truth, nearest track)", "Truth electrons", logy=True, atltext=atltext)
plotHistograms({f"|#eta| < {CEN_HI}": h["cen_dr_clu"], f"{FWD_LO} < |#eta| < {FWD_HI}": h["fwd_dr_clu"]},
               os.path.join(PLOT_DIR, "dr_cluster_by_region.png"), "#DeltaR(truth, nearest cluster)", "Truth electrons", logy=True, atltext=atltext)
plotHistograms({f"|#eta| < {CEN_HI}": h["cen_dr_pfo"], f"{FWD_LO} < |#eta| < {FWD_HI}": h["fwd_dr_pfo"]},
               os.path.join(PLOT_DIR, "dr_pfo_by_region.png"), "#DeltaR(truth, nearest PFO)", "Truth electrons", logy=True, atltext=atltext)
plotHistograms({f"|#eta| < {CEN_HI}": h["cen_Eratio"], f"{FWD_LO} < |#eta| < {FWD_HI}": h["fwd_Eratio"]},
               os.path.join(PLOT_DIR, "Eratio_by_region.png"), "#Sigma E^{clu}_{cone}/E^{truth}", "Truth electrons", atltext=atltext)
plotHistograms({f"|#eta| < {CEN_HI}": h["cen_nclu"], f"{FWD_LO} < |#eta| < {FWD_HI}": h["fwd_nclu"]},
               os.path.join(PLOT_DIR, "nclusters_by_region.png"), f"clusters within #DeltaR < {DR_CONE}", "Truth electrons", logy=True, atltext=atltext)
plotHistograms({f"|#eta| < {CEN_HI}": h["cen_clu_E"], f"{FWD_LO} < |#eta| < {FWD_HI}": h["fwd_clu_E"]},
               os.path.join(PLOT_DIR, "cluster_energy_by_region.png"), "Nearest cluster energy [GeV]", "Truth electrons", atltext=atltext)

eff_map = {}
for name, num in [("track", h["fwd_mcp_eta_trk"]), ("cluster", h["fwd_mcp_eta_clu"]), ("PFO", h["fwd_mcp_eta_pfo"])]:
    eff_map[name] = ROOT.TEfficiency(num, h["fwd_mcp_eta"]).CreateGraph()
plotEfficiencies(eff_map, os.path.join(PLOT_DIR, "forward_efficiency.png"), xlabel="|#eta|", ylabel="Efficiency")

print(f"\n\n--- {n_events} events, rates by |#eta| slice ---")
print(f"{'|eta|':<12} | {'N':<6} | {'track':<7} | {'cluster':<7} | {'PFO':<7} | {'e PFO':<7} | "
      f"{'chg PFO':<7} | {'trk+clu':<7} | {'<nClu>':<7} | {'<E/Etru>':<8}")
print("-"*106)
for lab in slice_labels:
    d = slices[lab]
    n = d["n"]
    if n == 0:
        print(f"{lab:<12} | {0:<6} |")
        continue
    mean_nclu = d["sum_nclu"]/n
    mean_er = d["sum_Eratio"]/d["n_Eratio"] if d["n_Eratio"] > 0 else 0
    print(f"{lab:<12} | {n:<6} | {d['trk']/n:<7.4f} | {d['clu']/n:<7.4f} | {d['pfo']/n:<7.4f} | "
          f"{d['pfo_el']/n:<7.4f} | {d['pfo_ch']/n:<7.4f} | {d['trkclu']/n:<7.4f} | {mean_nclu:<7.3f} | {mean_er:<8.4f}")

print(f"\n--- PFO composition, central (|#eta| < {CEN_HI}) vs forward ({FWD_LO} < |#eta| < {FWD_HI}) ---")
all_pdgs = sorted(set(list(pdg_counts["central"].keys()) + list(pdg_counts["forward"].keys())))
tot_cen = sum(pdg_counts["central"].values())
tot_fwd = sum(pdg_counts["forward"].values())
print(f"{'|PDG|':<8} | {'central':<9} | {'frac':<8} | {'forward':<9} | {'frac':<8}")
print("-"*52)
for pdg in all_pdgs:
    c = pdg_counts["central"].get(pdg, 0)
    fw = pdg_counts["forward"].get(pdg, 0)
    fc = c/tot_cen if tot_cen > 0 else 0
    ff = fw/tot_fwd if tot_fwd > 0 else 0
    print(f"{pdg:<8} | {c:<9} | {fc:<8.4f} | {fw:<9} | {ff:<8.4f}")

print(f"\nMean #DeltaR to nearest cluster, central : {h['cen_dr_clu'].GetMean():.4f}")
print(f"Mean #DeltaR to nearest cluster, forward : {h['fwd_dr_clu'].GetMean():.4f}")
print(f"Mean #DeltaR to nearest PFO,     central : {h['cen_dr_pfo'].GetMean():.4f}")
print(f"Mean #DeltaR to nearest PFO,     forward : {h['fwd_dr_pfo'].GetMean():.4f}")
print(f"Mean cone E / E truth,           central : {h['cen_Eratio'].GetMean():.4f}")
print(f"Mean cone E / E truth,           forward : {h['fwd_Eratio'].GetMean():.4f}")
print(f"\nHighest |#eta| with a track   : {h['trk_eta_all'].GetXaxis().GetBinCenter(h['trk_eta_all'].FindLastBinAbove(0)):.2f}")
print(f"Highest |#eta| with a cluster : {h['clu_eta_all'].GetXaxis().GetBinCenter(h['clu_eta_all'].FindLastBinAbove(0)):.2f}")
print(f"\nPlots saved to {PLOT_DIR}")
