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
 * @file BABAR_2004_S632399.cc
 *
 * @brief BABAR phi multiplicity and spectrum at Upsilon(4S)
 * 
 * @todo @test phi spectrum
 **/ 



namespace Rivet {
    /// @author Michael Große
    
    
    /**
     * @class BABAR_2004_S632399
     *
     * @brief BABAR phi spectrum at Upsilon(4S)
     *
     * @details This class is used to analyse the phi multiplicity at the Upsilon(4S) resonance.
     */
    class BABAR_2004_S632399 : public Analysis {
    public:
        
        /// The Constructor
        BABAR_2004_S632399()
        : Analysis("BABAR_2004_S632399"), _weightSum(0.)
        { }
        
    public:
        
        /**
         * @brief Books histograms and initialises counters and projections before the run
         *
         */
        void init() {
            _histMult   = bookHisto1D(1, 1, 1);
            _histdSigDp   = bookHisto1D(2, 1, 1);
            _phitotal=0;
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
            int const Y4S_id=300553;
            foreach (GenParticle* particle_possible_Y4S, Rivet::particles(e.genEvent())) {
                if(particle_possible_Y4S->pdg_id()!=Y4S_id) continue;
                GenVertex const * production_vertex = particle_possible_Y4S->production_vertex();
                bool passed = true;
                if (production_vertex) {
                    for (GenVertex::particles_in_const_iterator pp = production_vertex->particles_in_const_begin() ;pp != production_vertex->particles_in_const_end() ; ++pp) {
                        if ( particle_possible_Y4S->pdg_id() == (*pp)->pdg_id() ) {
                            passed = false;
                            break;
                        }
                    }
                }
                if(passed) {
                    upsilons.push_back(Particle(*particle_possible_Y4S));
                }
            }
            
            foreach (const Particle& p, upsilons) {
                vector<GenParticle *> phi;
                const int phiID = 333;
                findDecayProductsID(*p.genParticle(),phi,phiID);
                if (testBDecay(*p.genParticle())){
                    _weightSum += weight*2;
                }
                if((testBDecay(*p.genParticle()))&&(phi.size())){
                    LorentzTransform cms_boost(-p.momentum().boostVector());
                    for(unsigned int ix=0;ix<phi.size();++ix) {
                        double pcm = cms_boost.transform(FourMomentum(phi.at(ix)->momentum())).vector3().mod();
                        _histdSigDp->fill(pcm,weight);
                    }
                    _histMult->fill(0.,double(phi.size())*weight);
                    _phitotal+=phi.size();
                }
            }
        } // analyze
        
        /// Normalise histograms etc., after the run
        void finalize() {
            MSG_INFO("_phitotal " << _phitotal);
            scale(_histMult  , 1./_weightSum);
            scale(_histdSigDp  , 1./_weightSum);
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
        
        /// sum of the phi that were produced in B decays
        int _phitotal;
        
        
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
        * @param[in] upsilon The Upsilon(4S) to test if it decays into two B Mesons
        * @return This function returns \b true in case the Upsilons decay into two B-Mesons and \b false otherwise.
        *
        */
        bool testBDecay(GenParticle const & upsilon){
            int const particle_id = upsilon.pdg_id();
            int const Y4S_id=300553;
            if (particle_id!=Y4S_id){
                MSG_ERROR("Not an Upsilon(4S)");
                return false;
            }
            GenVertex const * end_vertex = upsilon.end_vertex();
            bool ret_is_B_Decay=true;
            for (GenVertex::particles_out_const_iterator daughter_iterator = end_vertex->particles_out_const_begin();\
                 daughter_iterator != end_vertex->particles_out_const_end(); ++daughter_iterator) {
                int const daughter_id = abs((*daughter_iterator)->pdg_id());
                int const Bplus_id = 521;
                int const Bzero_id = 511;
                if ((daughter_id!=Bplus_id)&&(daughter_id!=Bzero_id)){
                    ret_is_B_Decay=false;
                }
            }
            return ret_is_B_Decay;
        } //  testBDecay(...)
        
    }; // class
    
    
    /// The hook for the plugin system
    DECLARE_RIVET_PLUGIN(BABAR_2004_S632399);
    
} //namespace
