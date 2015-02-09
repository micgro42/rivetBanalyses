// -*- C++ -*-
#include <vector>
#include "Rivet/Analysis.hh"
//#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
//#include "Rivet/Projections/FinalState.hh"

/**
 * @file BABAR_gap.cc
 *
 * @brief looking at B -> X_c l nu
 *
 **/ 


namespace Rivet {
  /// @author Michael Große
  
  
  /**
   * @class BABAR_gap
   *
   * @brief looking at B -> X_c l nu
   *
   * @details detailed description
  */
  class BABAR_gap : public Analysis {
  public:

    /// The Constructor
    BABAR_gap()
    : Analysis("BABAR_gap"), _weightSum(0.)
    { }
    
  public:
    
    /**
     * @brief Book histograms and initialise counters and projections before the run
     *
     */
    void init() {
      
      
      _histDMult   = bookHisto1D(1, 1, 1);
      _histDstarMult   = bookHisto1D(2, 1, 1);
      _histD0starMult   = bookHisto1D(3, 1, 1);
      _histD1starMult   = bookHisto1D(4, 1, 1);
      _histD1Mult   = bookHisto1D(5, 1, 1);
      _histD2starMult   = bookHisto1D(6, 1, 1);
      _histDpiMult   = bookHisto1D(7, 1, 1);
      _histDstarpiMult   = bookHisto1D(8, 1, 1);
      _histXcMult   = bookHisto1D(9, 1, 1);
     
      
      
      
      // _h_XXXX = bookProfile1D(1, 1, 1);
      // _h_YYYY = bookHistogram1D(2, 1, 1);
      cout << "Initialising BABAR_gap done" << endl;
      
    }
    
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
              for (GenVertex::particles_in_const_iterator pp = pv->particles_in_const_begin() ;
                   pp != pv->particles_in_const_end() ; ++pp) {
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
        cout.precision(15);
        // find an upsilons
        foreach (const Particle& p, upsilons) {
          _weightSum += weight*2;
          // find B^- and /bar{B}^0
          vector<GenParticle *> B0;
          const int antiB0ID = -511;
          findDecayProducts(*p.genParticle(),B0,antiB0ID);
          _weightB0Sum += weight*B0.size();
          
          vector<GenParticle *> Bminus;
          const int BminusID = -521;
          findDecayProducts(*p.genParticle(),Bminus,BminusID);
          _weightBminusSum += weight*Bminus.size();
          
         
/// @todo check if D+ -> D+ + gamma => double counting
        const int D0ID = 421;
        const int D0barID = -421;
        const int DplusID = 411;
        const int DminusID = -411;
        const int DsplusID = 431;
        const int DsminusID = -431;
        const int LambdacplusID = 4122;
        const int LambdacminusID = -4122;
        if (Bminus.size()){
          //find $D^0$ from Bminus
          vector<GenParticle *> D0Bminus;
          findDecayProducts(*Bminus[0],D0Bminus,D0ID);
          _histD0BminusMult->fill(0.,double(D0Bminus.size())*weight);
          
          //find $\bar{D}^0$ from Bminus
          vector<GenParticle *> D0barBminus;
          findDecayProducts(*Bminus[0],D0barBminus,D0barID);
          _histD0barBminusMult->fill(0.,double(D0barBminus.size())*weight);
          
          //find $D^+$ from Bminus
          vector<GenParticle *> DpBminus;
          findDecayProducts(*Bminus[0],DpBminus,DplusID);
          _histDpBminusMult->fill(0.,double(DpBminus.size())*weight);
          
          //find $D^-$ from Bminus
          vector<GenParticle *> DmBminus;
          findDecayProducts(*Bminus[0],DmBminus,DminusID);
          _histDmBminusMult->fill(0.,double(DmBminus.size())*weight);
          
          //find $D_s^+$ from Bminus
          vector<GenParticle *> DspBminus;
          findDecayProducts(*Bminus[0],DspBminus,DsplusID);
          _histDspBminusMult->fill(0.,double(DspBminus.size())*weight);
          
          //find $D_s^-$ from Bminus
          vector<GenParticle *> DsmBminus;
          findDecayProducts(*Bminus[0],DsmBminus,DsminusID);
          _histDsmBminusMult->fill(0.,double(DsmBminus.size())*weight);
          
          //find $\Lambda_c^+$ from Bminus
          vector<GenParticle *> LambdacpBminus;
          findDecayProducts(*Bminus[0],LambdacpBminus,LambdacplusID);
          _histLambdacpBminusMult->fill(0.,double(LambdacpBminus.size())*weight);
          
          //find $\Lambda_c^-$ from Bminus
          vector<GenParticle *> LambdacmBminus;
          findDecayProducts(*Bminus[0],LambdacmBminus,LambdacminusID);
          _histLambdacmBminusMult->fill(0.,double(LambdacmBminus.size())*weight);
        }
        
        if (B0.size()){
          //find $D^0$ from B0
          vector<GenParticle *> D0B0;
          findDecayProducts(*B0[0],D0B0,D0ID);
          _histD0B0Mult->fill(0.,double(D0B0.size())*weight);
          
          //find $\bar{D}^0$ from B0
          vector<GenParticle *> D0barB0;
          findDecayProducts(*B0[0],D0barB0,D0barID);
          _histD0barB0Mult->fill(0.,double(D0barB0.size())*weight);
          
          //find $D^+$ from B0
          vector<GenParticle *> DpB0;
          findDecayProducts(*B0[0],DpB0,DplusID);
          _histDpB0Mult->fill(0.,double(DpB0.size())*weight);
          
          //find $D^-$ from B0
          vector<GenParticle *> DmB0;
          findDecayProducts(*B0[0],DmB0,DminusID);
          _histDmB0Mult->fill(0.,double(DmB0.size())*weight);
          
          //find $D_s^+$ from B0
          vector<GenParticle *> DspB0;
          findDecayProducts(*B0[0],DspB0,DsplusID);
          _histDspB0Mult->fill(0.,double(DspB0.size())*weight);
          
          //find $D_s^-$ from B0
          vector<GenParticle *> DsmB0;
          findDecayProducts(*B0[0],DsmB0,DsminusID);
          _histDsmB0Mult->fill(0.,double(DsmB0.size())*weight);
          
          //find $\Lambda_c^+$ from B0
          vector<GenParticle *> LambdacpB0;
          findDecayProducts(*B0[0],LambdacpB0,LambdacplusID);
          _histLambdacpB0Mult->fill(0.,double(LambdacpB0.size())*weight);
          
          //find $\Lambda_c^-$ from B0
          vector<GenParticle *> LambdacmB0;
          findDecayProducts(*B0[0],LambdacmB0,LambdacminusID);
          _histLambdacmB0Mult->fill(0.,double(LambdacmB0.size())*weight);
        }
    }
  } // analyze
  
  
  
  /// Normalise histograms etc., after the run
  void finalize() {
    
    
    scale(_histD0BminusMult  , 1./_weightBminusSum);
    scale(_histD0barBminusMult  , 1./_weightBminusSum);
    scale(_histDpBminusMult  , 1./_weightBminusSum);
    scale(_histDmBminusMult  , 1./_weightBminusSum);
    scale(_histDspBminusMult  , 1./_weightBminusSum);
    scale(_histDsmBminusMult  , 1./_weightBminusSum);
    scale(_histLambdacpBminusMult  , 1./_weightBminusSum);
    scale(_histLambdacmBminusMult  , 1./_weightBminusSum);
    scale(_histD0B0Mult  , 1./_weightB0Sum);
    scale(_histD0barB0Mult  , 1./_weightB0Sum);
    scale(_histDpB0Mult  , 1./_weightB0Sum);
    scale(_histDmB0Mult  , 1./_weightB0Sum);
    scale(_histDspB0Mult  , 1./_weightB0Sum);
    scale(_histDsmB0Mult  , 1./_weightB0Sum);
    scale(_histLambdacpB0Mult  , 1./_weightB0Sum);
    scale(_histLambdacmB0Mult  , 1./_weightB0Sum);
    
    
  }
  
  //@}
  
  
private:
  
  
  
  /// sum of weights
  double _weightSum;
  
  
  /**
  * @name Histograms
  * 
   * @see doi:10.1103/PhysRevD.75.072002 Table III 
  */
  //@{
  
  
  
  /// Branching-ratio \f$ B \rightarrow D + l \overline{\nu} \f$
  Histo1DPtr _histDMult;
  
  /// Branching-ratio \f$ B \rightarrow D^{*} + l \overline{\nu} \f$
  Histo1DPtr _histDstarMult;
  
  /// Branching-ratio \f$ B \rightarrow D^{*}_{0} + l \overline{\nu} \f$
  Histo1DPtr _histD0starMult;
  
  /// Branching-ratio \f$ B \rightarrow D^{*}_{1} + l \overline{\nu} \f$
  Histo1DPtr _histD1starMult;
  
  /// Branching-ratio \f$ B \rightarrow D_{1} + l \overline{\nu} \f$
  Histo1DPtr _histD1Mult;
  
  /// Branching-ratio \f$ B \rightarrow D^{*}_{2} + l \overline{\nu} \f$
  Histo1DPtr _histD2starMult;
  
  /// Branching-ratio \f$ B \rightarrow D + \pi + l \overline{\nu} \f$
  Histo1DPtr _histDpiMult;
  
  /// Branching-ratio \f$ B \rightarrow D^{*} + \pi + l \overline{\nu} \f$
  Histo1DPtr _histDstarpiMult;
  
  /// Branching-ratio \f$ B \rightarrow X_{c} + l \overline{\nu} \f$
  Histo1DPtr _histXcMult;
  
  
  
 
  //@}
  
  /**
   * @brief function to find decayproducts by ID
   * 
   * @param[in] p the particle whose decay-tree will be searched
   * @param[in] selectionid the id of the particle(decay product) according to evt.pdl that is searched for
   * @param[out] pvec the decay products found will be written to this vector
   *
   *
   */ 
  void findDecayProducts(const GenParticle & p,
                         vector<GenParticle *> & pvec,int selectionid) {
    const GenVertex* dv = p.end_vertex();
    for (GenVertex::particles_out_const_iterator
      pp = dv->particles_out_const_begin();
    pp != dv->particles_out_const_end(); ++pp) {
      int id = (*pp)->pdg_id();
      if(id==selectionid) { //rivet.hepforge.org/svn/trunk/include/Rivet/ParticleName.hh
        pvec.push_back(*pp);
      }
      else if((*pp)->end_vertex())
        findDecayProducts(**pp,pvec,selectionid);
    }
  }
                         
  ///blob
   void printDaughters(vector<int> DaughterIDs){
       for (std::vector<int>::iterator it = DaughterIDs.begin() ; it != DaughterIDs.end(); ++it) {
                             cout << *it << " ";
                           }
                           cout << endl;
                           
                         }


};



/// The hook for the plugin system
DECLARE_RIVET_PLUGIN(BABAR_gap);

}
