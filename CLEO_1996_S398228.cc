// -*- C++ -*-
#include <iostream>
#include <vector>
#include "Rivet/Analysis.hh"
//#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Tools/Logging.hh"

/**
 * @file CLEO_1996_S398228.cc
 *
 * @brief CLEO \f$ \eta \f$ multiplicity and spectrum at Upsilon(4S)
 * 
 **/ 



namespace Rivet {
    /// @author Michael Große
    
    
    /**
     * @class CLEO_1996_S398228
     *
     * @brief BABAR \f$ \eta \f$ spectrum at Upsilon(4S)
     *
     * @details This class is used to analyse the \f$ \eta \f$ multiplicity and momentum spectrum at the Upsilon(4S) resonance.
     */
    class CLEO_1996_S398228 : public Analysis {
    public:
        
        /// The Constructor
        CLEO_1996_S398228()
            : Analysis("CLEO_1996_S398228"), _weightSum(0.)
        { }
        
    public:
        
        /**
         * @brief Books histograms and initialises counters and projections before the run
         *
         */
        void init() {
            _histMult   = bookHisto1D(1, 1, 1);
            _histdSigDp   = bookHisto1D(2, 1, 1);
            _etatotal=0;
        } // init
        
        /** 
         * @brief Perform the per-event analysis
         *
         * @param[in] e event created by the MC Generator
         *
         */
        void analyze(const Event& e) {
            const double weight = e.weight();
            
            // find the upsilons
            ParticleVector upsilons;
            foreach (GenParticle* p, Rivet::particles(e.genEvent())) {
                if(p->pdg_id()!=300553) continue;
                const GenVertex* pv = p->production_vertex();
                bool passed = true;
                if (pv) {
                    for (GenVertex::particles_in_const_iterator pp = pv->particles_in_const_begin() ;pp != pv->particles_in_const_end() ; ++pp) {
                        if ( p->pdg_id() == (*pp)->pdg_id() ) {
                            passed = false;
                            break;
                        }
                    }
                }
                if(passed) {
                    upsilons.push_back(Particle(*p));
                }
            }
            
            foreach (const Particle& p, upsilons) {
                
                vector<GenParticle *> eta;
                const int etaID = 221;
                findDecayProductsID(*p.genParticle(),eta,etaID);
                if (testBDecay(*p.genParticle())){
                    _weightSum += weight*2;
                }
                if((testBDecay(*p.genParticle()))&&(eta.size())){
                    LorentzTransform cms_boost(-p.momentum().boostVector());
                    for(unsigned int ix=0;ix<eta.size();++ix) {
                        double pcm = cms_boost.transform(FourMomentum(eta.at(ix)->momentum())).vector3().mod();
                        _histdSigDp->fill(pcm,weight);
                    }
                    _histMult->fill(0.,(double)eta.size()*weight);
                    _etatotal+=eta.size()*weight;
                }
            }
        } // analyze
        
        /// Normalise histograms etc., after the run
        void finalize() {
            MSG_INFO("_etatotal " << _etatotal);
            scale(_histMult  , 1./_weightSum);
            scale(_histdSigDp  , 1./(10*_etatotal));
            MSG_INFO("_weightSum, i.e. number of B mesons = " << _weightSum);
            
        } // finalize
        
        
        
        
    private:
        
        
        
        /// @name Histograms
        //@{
        Histo1DPtr _histMult;
        Histo1DPtr _histdSigDp;
        //@}
        
        /// sum of weights, i. e. how may B Mesons were produced
        double _weightSum;
        
        /// sum of the eta that were produced in B decays
        int _etatotal;
        
        
        /**
         * @brief function to find decayproducts by ID
         * 
         * @param[in] p the particle whose decay-tree will be searched
         * @param[in] selectionid the id of the particle(decay product) according to evt.pdl that is searched for
         * @param[out] pvec the decay products found will be written to this vector
         *
         *
        */
        void findDecayProductsID(const GenParticle & p, vector<GenParticle *> & pvec,int selectionid) {
            const GenVertex* dv = p.end_vertex();
            for (GenVertex::particles_out_const_iterator pp = dv->particles_out_const_begin(); pp != dv->particles_out_const_end(); ++pp) {
                int id = abs((*pp)->pdg_id());
                if(id==selectionid) { 
                    pvec.push_back(*pp);
                }else if((*pp)->end_vertex()){
                    findDecayProductsID(**pp,pvec,selectionid);
                }
            }
        }
        
        
        /**
        * 
        * @brief Test if the Upsilon(4S) Decays into B Mesons
        * @param[in] p The Upsilon(4S) to test if it decays into two B Mesons
        * @return This function returns \b true in case the Upsilons decay into two B-Mesons and \b false otherwise.
        *
        */
        bool testBDecay(const GenParticle & p){
            int mid = p.pdg_id();
            if (mid!=300553){
                MSG_ERROR("Not a Upsilon(4S)");
                return false;
            }
            const GenVertex* dv = p.end_vertex();
            bool ret_is_B_Decay=true;
            for (GenVertex::particles_out_const_iterator pp = dv->particles_out_const_begin(); pp != dv->particles_out_const_end(); ++pp) {
                int id = abs((*pp)->pdg_id());
                const int BplusID = 521;
                const int BzeroID = 511;
                if ((id!=BplusID)&&(id!=BzeroID)){
                    ret_is_B_Decay=false;
                }
            }
            return ret_is_B_Decay;
        } //  testBDecay(...)
        
    }; // class
    
    
    /// The hook for the plugin system
    DECLARE_RIVET_PLUGIN(CLEO_1996_S398228);
    
} //namespace
