import math
import glob
import ROOT
import os
from podio import root_io
exec(open("./plotHelper.py").read())
ROOT.gROOT.SetBatch()
PLOT_DIR = "/scratch/jwatts/mucol/mucolstudies/plots2026/diagnostics"
os.makedirs(PLOT_DIR, exist_ok=True)
samples = sorted(glob.glob("/scratch/jwatts/mucol/v3/reco/electronGun_0_50/*_reco_*.edm4hep.root"))
#samples = sorted(glob.glob("/scratch/jwatts/mucol/v3/reco/electronGun_50_250/*_reco_*.edm4hep.root"))
#samples = sorted(glob.glob("/scratch/jwatts/mucol/v3/reco/electronGun_250_1000/*_reco_*.edm4hep.root"))
#samples = sorted(glob.glob("/scratch/jwatts/mucol/v3/reco/electronGun_1000_5000/*_reco_*.edm4hep.root"))
files = {"electronGun_0_50": samples}

B_FIELD = 5
M_ELECTRON = 0.000511
DR_CONE = 0.2
PT_MIN = 2.0            # GeV, track threshold for the per event table
E_MIN = 2.0             # GeV, cluster and PFO threshold for the per event table
ETA_MAX = 2.4
N_PRINT_EVENTS = 20     # how many events to print one line each for

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

# counters for the summary table, one increment per event
counts = {}
def countIf(label, condition):
    if label not in counts: counts[label] = 0
    if condition: counts[label] += 1
    return

# running sums for the mean-per-event table
sums = {}
def addTo(label, value):
    if label not in sums: sums[label] = 0
    sums[label] += value
    return

h = {}
# per event multiplicities
h["n_mcp"] = ROOT.TH1F("n_mcp", "", 21, -0.5, 20.5)
h["n_trk"] = ROOT.TH1F("n_trk", "", 21, -0.5, 20.5)
h["n_clu"] = ROOT.TH1F("n_clu", "", 21, -0.5, 20.5)
h["n_pfo"] = ROOT.TH1F("n_pfo", "", 21, -0.5, 20.5)
h["n_pfo_el"] = ROOT.TH1F("n_pfo_el", "", 21, -0.5, 20.5)
# truth kinematics
h["mcp_pt"]  = ROOT.TH1F("mcp_pt",  "", 60, 0, 60)
h["mcp_eta"] = ROOT.TH1F("mcp_eta", "", 80, -4, 4)
h["mcp_phi"] = ROOT.TH1F("mcp_phi", "", 70, -3.5, 3.5)
h["mcp_E"]   = ROOT.TH1F("mcp_E",   "", 60, 0, 60)
# reco kinematics, every object in the collection, no matching
h["trk_pt"]  = ROOT.TH1F("trk_pt",  "", 60, 0, 60)
h["trk_eta"] = ROOT.TH1F("trk_eta", "", 80, -4, 4)
h["trk_phi"] = ROOT.TH1F("trk_phi", "", 70, -3.5, 3.5)
h["clu_E"]   = ROOT.TH1F("clu_E",   "", 60, 0, 60)
h["clu_eta"] = ROOT.TH1F("clu_eta", "", 80, -4, 4)
h["clu_phi"] = ROOT.TH1F("clu_phi", "", 70, -3.5, 3.5)
h["pfo_E"]   = ROOT.TH1F("pfo_E",   "", 60, 0, 60)
h["pfo_eta"] = ROOT.TH1F("pfo_eta", "", 80, -4, 4)
h["pfo_phi"] = ROOT.TH1F("pfo_phi", "", 70, -3.5, 3.5)
# track quality
h["trk_chi2ndf"] = ROOT.TH1F("trk_chi2ndf", "", 50, 0, 20)
h["trk_ndf"]     = ROOT.TH1F("trk_ndf",     "", 41, -0.5, 40.5)
# how far away is the nearest reco object from the truth electron
h["dr_trk"] = ROOT.TH1F("dr_trk", "", 100, 0, 2)
h["dr_clu"] = ROOT.TH1F("dr_clu", "", 100, 0, 2)
h["dr_pfo"] = ROOT.TH1F("dr_pfo", "", 100, 0, 2)
h["dr_pfo_any"] = ROOT.TH1F("dr_pfo_any", "", 100, 0, 2)
# how many objects land inside the cone
h["n_trk_cone"] = ROOT.TH1F("n_trk_cone", "", 11, -0.5, 10.5)
h["n_clu_cone"] = ROOT.TH1F("n_clu_cone", "", 11, -0.5, 10.5)
h["n_pfo_cone"] = ROOT.TH1F("n_pfo_cone", "", 11, -0.5, 10.5)
# response of the nearest matched object
h["trk_dpt"]   = ROOT.TH1F("trk_dpt",   "", 100, -1, 1)
h["clu_Eratio"] = ROOT.TH1F("clu_Eratio", "", 100, 0, 2)
h["pfo_Eratio"] = ROOT.TH1F("pfo_Eratio", "", 100, 0, 2)
# truth eta of electrons that do or do not get a track, to see where it dies
h["mcp_eta_all"]    = ROOT.TH1F("mcp_eta_all",    "", 50, -2.5, 2.5)
h["mcp_eta_notrk"]  = ROOT.TH1F("mcp_eta_notrk",  "", 50, -2.5, 2.5)
h["mcp_eta_noclu"]  = ROOT.TH1F("mcp_eta_noclu",  "", 50, -2.5, 2.5)
h["mcp_pt_all"]     = ROOT.TH1F("mcp_pt_all",     "", 60, 0, 60)
h["mcp_pt_notrk"]   = ROOT.TH1F("mcp_pt_notrk",   "", 60, 0, 60)

pfo_pdg_counts = {}
n_events = 0
resolved = False

print(f"\n--- First {N_PRINT_EVENTS} events ---")
print(f"{'Event':<6} | {'nTruth':<7} | {'nTrk':<5} | {'pT>2':<5} | {'nClu':<5} | {'E>2':<5} | {'nPFO':<5} | {'nElPFO':<6}")

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
                h["mcp_pt"].Fill(tlv.Perp())
                h["mcp_eta"].Fill(tlv.Eta())
                h["mcp_phi"].Fill(tlv.Phi())
                h["mcp_E"].Fill(tlv.E())
            mcp_central = [t for t in mcp_electrons if abs(t.Eta()) < ETA_MAX]

            trk_tlvs = []
            for t in trks:
                tlv = getTrackTLV_edm4hep(t)
                trk_tlvs.append(tlv)
                h["trk_pt"].Fill(tlv.Perp())
                h["trk_eta"].Fill(tlv.Eta())
                h["trk_phi"].Fill(tlv.Phi())
                ndf = t.getNdf()
                h["trk_ndf"].Fill(ndf)
                if ndf > 0: h["trk_chi2ndf"].Fill(t.getChi2()/ndf)

            clu_tlvs = []
            for c in clusters:
                tlv = getClusterTLV_edm4hep(c)
                clu_tlvs.append(tlv)
                h["clu_E"].Fill(tlv.E())
                h["clu_eta"].Fill(tlv.Eta())
                h["clu_phi"].Fill(tlv.Phi())

            pfo_tlvs_all = []
            pfo_tlvs_el = []
            pfo_tlvs_charged = []
            for p in pfos:
                tlv = getRecoTLV(p)
                pdg = abs(getPDG(p))
                pfo_tlvs_all.append(tlv)
                pfo_pdg_counts[pdg] = pfo_pdg_counts.get(pdg, 0) + 1
                if pdg == 11: pfo_tlvs_el.append(tlv)
                if pdg in (11, 211, 13): pfo_tlvs_charged.append(tlv)
                h["pfo_E"].Fill(tlv.E())
                h["pfo_eta"].Fill(tlv.Eta())
                h["pfo_phi"].Fill(tlv.Phi())

            trk_above = [t for t in trk_tlvs if t.Perp() > PT_MIN]
            clu_above = [c for c in clu_tlvs if c.E() > E_MIN]
            pfo_el_above = [p for p in pfo_tlvs_el if p.E() > E_MIN]

            h["n_mcp"].Fill(len(mcp_electrons))
            h["n_trk"].Fill(len(trk_tlvs))
            h["n_clu"].Fill(len(clu_tlvs))
            h["n_pfo"].Fill(len(pfo_tlvs_all))
            h["n_pfo_el"].Fill(len(pfo_tlvs_el))

            if n_events <= N_PRINT_EVENTS:
                print(f"{n_events:<6} | {len(mcp_electrons):<7} | {len(trk_tlvs):<5} | {len(trk_above):<5} | "
                      f"{len(clu_tlvs):<5} | {len(clu_above):<5} | {len(pfo_tlvs_all):<5} | {len(pfo_tlvs_el):<6}")

            # per event flags, matched versions get set inside the truth loop below
            cone = {"trk": False, "trk_pt": False, "clu": False, "clu_E": False,
                    "pfo": False, "pfo_el": False, "pfo_charged": False}

            for mcp_el in mcp_electrons:
                m_eta = mcp_el.Eta()
                m_pt = mcp_el.Perp()
                h["mcp_eta_all"].Fill(m_eta)
                h["mcp_pt_all"].Fill(m_pt)

                drs_trk = sorted([(mcp_el.DeltaR(t), t) for t in trk_tlvs], key=lambda x: x[0])
                drs_clu = sorted([(mcp_el.DeltaR(c), c) for c in clu_tlvs], key=lambda x: x[0])
                drs_pfo = sorted([(mcp_el.DeltaR(p), p) for p in pfo_tlvs_el], key=lambda x: x[0])
                drs_pfo_any = sorted([(mcp_el.DeltaR(p), p) for p in pfo_tlvs_all], key=lambda x: x[0])
                drs_pfo_ch = sorted([(mcp_el.DeltaR(p), p) for p in pfo_tlvs_charged], key=lambda x: x[0])

                if len(drs_trk) > 0: h["dr_trk"].Fill(drs_trk[0][0])
                if len(drs_clu) > 0: h["dr_clu"].Fill(drs_clu[0][0])
                if len(drs_pfo) > 0: h["dr_pfo"].Fill(drs_pfo[0][0])
                if len(drs_pfo_any) > 0: h["dr_pfo_any"].Fill(drs_pfo_any[0][0])

                n_trk_cone = len([1 for dr, t in drs_trk if dr < DR_CONE])
                n_clu_cone = len([1 for dr, c in drs_clu if dr < DR_CONE])
                n_pfo_cone = len([1 for dr, p in drs_pfo if dr < DR_CONE])
                h["n_trk_cone"].Fill(n_trk_cone)
                h["n_clu_cone"].Fill(n_clu_cone)
                h["n_pfo_cone"].Fill(n_pfo_cone)

                if n_trk_cone > 0: cone["trk"] = True
                if n_clu_cone > 0: cone["clu"] = True
                if n_pfo_cone > 0: cone["pfo_el"] = True
                if len([1 for dr, t in drs_trk if dr < DR_CONE and t.Perp() > PT_MIN]) > 0: cone["trk_pt"] = True
                if len([1 for dr, c in drs_clu if dr < DR_CONE and c.E() > E_MIN]) > 0: cone["clu_E"] = True
                if len([1 for dr, p in drs_pfo_any if dr < DR_CONE]) > 0: cone["pfo"] = True
                if len([1 for dr, p in drs_pfo_ch if dr < DR_CONE]) > 0: cone["pfo_charged"] = True

                if n_trk_cone > 0 and m_pt > 0:
                    h["trk_dpt"].Fill((drs_trk[0][1].Perp() - m_pt)/m_pt)
                else:
                    h["mcp_eta_notrk"].Fill(m_eta)
                    h["mcp_pt_notrk"].Fill(m_pt)
                if n_clu_cone > 0 and mcp_el.E() > 0:
                    h["clu_Eratio"].Fill(drs_clu[0][1].E()/mcp_el.E())
                else:
                    h["mcp_eta_noclu"].Fill(m_eta)
                if n_pfo_cone > 0 and mcp_el.E() > 0:
                    h["pfo_Eratio"].Fill(drs_pfo[0][1].E()/mcp_el.E())

            # anywhere in the event, no matching required
            countIf("truth electron in event",                 len(mcp_electrons) > 0)
            countIf(f"truth electron, |#eta| < {ETA_MAX}",      len(mcp_central) > 0)
            countIf("any track",                               len(trk_tlvs) > 0)
            countIf(f"any track, pT > {PT_MIN:.0f} GeV",        len(trk_above) > 0)
            countIf("any cluster",                             len(clu_tlvs) > 0)
            countIf(f"any cluster, E > {E_MIN:.0f} GeV",        len(clu_above) > 0)
            countIf("any PFO",                                 len(pfo_tlvs_all) > 0)
            countIf("any electron PFO",                        len(pfo_tlvs_el) > 0)
            countIf(f"any electron PFO, E > {E_MIN:.0f} GeV",   len(pfo_el_above) > 0)
            countIf("any charged PFO (11, 13, 211)",           len(pfo_tlvs_charged) > 0)
            # matched to the truth electron within the cone
            countIf(f"track within #DeltaR < {DR_CONE}",        cone["trk"])
            countIf(f"track in cone, pT > {PT_MIN:.0f} GeV",    cone["trk_pt"])
            countIf(f"cluster within #DeltaR < {DR_CONE}",      cone["clu"])
            countIf(f"cluster in cone, E > {E_MIN:.0f} GeV",    cone["clu_E"])
            countIf("PFO in cone (any type)",                  cone["pfo"])
            countIf("electron PFO in cone",                    cone["pfo_el"])
            countIf("charged PFO in cone",                     cone["pfo_charged"])
            countIf("track AND cluster in cone",               cone["trk"] and cone["clu"])

            addTo("truth electrons",   len(mcp_electrons))
            addTo("tracks",            len(trk_tlvs))
            addTo(f"tracks pT > {PT_MIN:.0f} GeV", len(trk_above))
            addTo("clusters",          len(clu_tlvs))
            addTo(f"clusters E > {E_MIN:.0f} GeV", len(clu_above))
            addTo("PFOs",              len(pfo_tlvs_all))
            addTo("electron PFOs",     len(pfo_tlvs_el))
            addTo("charged PFOs",      len(pfo_tlvs_charged))

atltext = ["Muon Collider", "Simulation, no BIB"]

plotHistograms({"truth electrons": h["n_mcp"], "tracks": h["n_trk"], "clusters": h["n_clu"], "PFOs": h["n_pfo"], "electron PFOs": h["n_pfo_el"]},
               os.path.join(PLOT_DIR, "multiplicity.png"), "number per event", "Events", logy=True, atltext=atltext)
plotHistograms({"truth electron": h["mcp_pt"], "track": h["trk_pt"]},
               os.path.join(PLOT_DIR, "pt_truth_vs_track.png"), "p_{T} [GeV]", "Objects", atltext=atltext)
plotHistograms({"truth electron": h["mcp_E"], "cluster": h["clu_E"], "PFO": h["pfo_E"]},
               os.path.join(PLOT_DIR, "energy_truth_vs_reco.png"), "Energy [GeV]", "Objects", atltext=atltext)
plotHistograms({"truth electron": h["mcp_eta"], "track": h["trk_eta"], "cluster": h["clu_eta"], "PFO": h["pfo_eta"]},
               os.path.join(PLOT_DIR, "eta_truth_vs_reco.png"), "#eta", "Objects", atltext=atltext)
plotHistograms({"truth electron": h["mcp_phi"], "track": h["trk_phi"], "cluster": h["clu_phi"], "PFO": h["pfo_phi"]},
               os.path.join(PLOT_DIR, "phi_truth_vs_reco.png"), "#phi", "Objects", atltext=atltext)
plotHistograms({"#chi^{2}/ndf": h["trk_chi2ndf"]},
               os.path.join(PLOT_DIR, "track_chi2ndf.png"), "#chi^{2}/ndf", "Tracks", atltext=atltext)
plotHistograms({"ndf": h["trk_ndf"]},
               os.path.join(PLOT_DIR, "track_ndf.png"), "ndf", "Tracks", atltext=atltext)
plotHistograms({"nearest track": h["dr_trk"], "nearest cluster": h["dr_clu"],
                "nearest PFO (PDG=11)": h["dr_pfo"], "nearest PFO (any)": h["dr_pfo_any"]},
               os.path.join(PLOT_DIR, "dr_nearest.png"), "#DeltaR(truth, reco)", "Truth electrons", logy=True, atltext=atltext)
plotHistograms({"tracks": h["n_trk_cone"], "clusters": h["n_clu_cone"], "PFOs (PDG=11)": h["n_pfo_cone"]},
               os.path.join(PLOT_DIR, "n_in_cone.png"), f"number within #DeltaR < {DR_CONE}", "Truth electrons", logy=True, atltext=atltext)
plotHistograms({"(p_{T}^{trk} - p_{T}^{truth})/p_{T}^{truth}": h["trk_dpt"]},
               os.path.join(PLOT_DIR, "track_pt_response.png"), "(p_{T}^{trk} - p_{T}^{truth})/p_{T}^{truth}", "Matches", atltext=atltext)
plotHistograms({"cluster": h["clu_Eratio"], "PFO": h["pfo_Eratio"]},
               os.path.join(PLOT_DIR, "energy_response.png"), "E^{reco}/E^{truth}", "Matches", atltext=atltext)
plotHistograms({"all truth electrons": h["mcp_eta_all"], "no track in cone": h["mcp_eta_notrk"], "no cluster in cone": h["mcp_eta_noclu"]},
               os.path.join(PLOT_DIR, "eta_unmatched.png"), "#eta", "Truth electrons", atltext=atltext)
plotHistograms({"all truth electrons": h["mcp_pt_all"], "no track in cone": h["mcp_pt_notrk"]},
               os.path.join(PLOT_DIR, "pt_unmatched.png"), "p_{T} [GeV]", "Truth electrons", atltext=atltext)

print(f"\n--- {n_events} events total ---")
print(f"collections: mcp={coll_names['mcp']}, trk={coll_names['trk']}, cluster={coll_names['cluster']}, pfo={coll_names['pfo']}\n")

print(f"{'Objects per event':<32} | {'Total':<10} | {'Mean':<8}")
print("-"*56)
for label in sums:
    print(f"{label:<32} | {sums[label]:<10} | {sums[label]/n_events:<8.3f}")

print(f"\n{'Event rate':<38} | {'Events':<8} | {'Rate':<8} | {'Error':<8}")
print("-"*70)
for label in counts:
    n = counts[label]
    rate = n/n_events if n_events > 0 else 0
    err = math.sqrt(rate*(1 - rate)/n_events) if n_events > 0 else 0
    print(f"{label:<38} | {n:<8} | {rate:<8.4f} | {err:<8.4f}")

print(f"\n{'PFO |PDG|':<15} | {'Count':<10} | {'Fraction':<10} | {'Per event':<10}")
print("-"*54)
total_pfos = sum(pfo_pdg_counts.values())
for pdg in sorted(pfo_pdg_counts, key=lambda k: -pfo_pdg_counts[k]):
    frac = pfo_pdg_counts[pdg]/total_pfos if total_pfos > 0 else 0
    print(f"{pdg:<15} | {pfo_pdg_counts[pdg]:<10} | {frac:<10.4f} | {pfo_pdg_counts[pdg]/n_events:<10.3f}")

print(f"\nMean #DeltaR to nearest track   : {h['dr_trk'].GetMean():.4f}")
print(f"Mean #DeltaR to nearest cluster : {h['dr_clu'].GetMean():.4f}")
print(f"Mean #DeltaR to nearest PFO     : {h['dr_pfo_any'].GetMean():.4f}")
print(f"Mean track p_T response         : {h['trk_dpt'].GetMean():.4f}")
print(f"Mean cluster E ratio            : {h['clu_Eratio'].GetMean():.4f}")
print(f"\nPlots saved to {PLOT_DIR}")
