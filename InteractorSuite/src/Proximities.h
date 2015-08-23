//
//  CallInteractions.h
//  InteractorSuite
//
//  Created by Pelin Sahlen on 13/08/2015.
//  Copyright (c) 2015 Pelin Sahlen. All rights reserved.
//
class ProximityClass{
    friend class ProcessBAM;
    friend class PromoterClass;
void AnnotateDistalInteractor(ProbeSet&, int, string, string, int, int, int);
void AnnotateFeatFeatInteraction(ProbeSet&, int, int, int);
void PopulateInteractions(boost::unordered::unordered_map<int, int* >&, int, int);
};

void ProximityClass::AnnotateDistalInteractor(ProbeSet& probes, int probe_index, string anchored_chr, string interactor_chr, int anchored_resite, int interactor_resite, int ExperimentNo){
    const int offset = 0;
    // Intra chromosomal interaction
    if(anchored_chr.compare(interactor_chr) == 0){
    //    cout << "dist " << abs(anchored_resite - interactor_resite) << "  " << probe_index << endl;
        if( (abs(anchored_resite - interactor_resite)) > MinimumJunctionDistance){ // at least minjunctdist away
            PopulateInteractions(probes.Probes[probe_index].feature.Signals.signals, interactor_resite, ExperimentNo);
        }
    }
    else{  // inter chromosomal interaction
        bool chrfound = 0;
        vector < SignalStruct_CTX >::iterator itx;
        for(itx = probes.Probes[probe_index].feature.Signals_CTX.begin() ; itx < probes.Probes[probe_index].feature.Signals_CTX.end();++itx){
            if (interactor_chr.compare(itx->maptochrname) == 0){
                if(itx->signal.find(interactor_resite) == itx->signal.end()){
                    itx->signal[interactor_resite] = new int[NOFEXPERIMENTS + offset]; // add a new entry
                    for(int z = 0; z < NOFEXPERIMENTS; ++z)
                        itx->signal[interactor_resite][z + offset] = 0;
                    
                    itx->signal[interactor_resite][ExperimentNo + offset] = 1;
                }
                else{ // if inserted before
                    itx->signal[interactor_resite][ExperimentNo + offset] += 1;
                }
                chrfound = 1;
                break;
            }
        }
        if(!chrfound){
            probes.Probes[probe_index].feature.Signals_CTX.push_back(SignalStruct_CTX());
            probes.Probes[probe_index].feature.Signals_CTX.back().maptochrname.append(interactor_chr);
            probes.Probes[probe_index].feature.Signals_CTX.back().signal[interactor_resite] = new int[NOFEXPERIMENTS + offset];
            for(int z = 0; z < NOFEXPERIMENTS; ++z)
                probes.Probes[probe_index].feature.Signals_CTX.back().signal[interactor_resite][z + offset] = 0;
            probes.Probes[probe_index].feature.Signals_CTX.back().signal[interactor_resite][ExperimentNo + offset] = 1;
        }
    }
//    cout << anchored_resite << "  " << interactor_resite <<  "  distal interactor annotated " << endl;
}

void ProximityClass::AnnotateFeatFeatInteraction(ProbeSet& probes, int probe_index1, int probe_index2, int ExperimentNo){
    bool foundbefore = 0;
    //Check if the second promoter is close
    vector < FeattoFeatSignalStruct >::iterator it;
    for (it = probes.Probes[probe_index1].feature.Inter_feature_ints.begin(); it < probes.Probes[probe_index1].feature.Inter_feature_ints.end(); ++it){ // Check if the interaction with that promoter is seen before
        if (it->feat_index == probe_index2){
            it->normsignal[ExperimentNo] += 1.0;
            foundbefore = 1;
            break;
        }
    }
    if(!foundbefore){ // Create a new entry for that promoter
        probes.Probes[probe_index1].feature.Inter_feature_ints.push_back(FeattoFeatSignalStruct());
        probes.Probes[probe_index1].feature.Inter_feature_ints.back().feat_index = probe_index2;
        probes.Probes[probe_index1].feature.Inter_feature_ints.back().normsignal[ExperimentNo] += 1.0;
    }
}

void ProximityClass::PopulateInteractions(boost::unordered::unordered_map<int, int* >& signals, int interactor_resite, int ExperimentNo){
    const int offset = 0; //[0] = interactor pos
    if(signals.find(interactor_resite) == signals.end()){
        signals[interactor_resite] = new int[NOFEXPERIMENTS + offset]; // add a new entry // [0] = interactor pos, [1] = promoter_REsite
        for(int z = 0; z < (NOFEXPERIMENTS + offset); ++z)
            signals[interactor_resite][z + offset] = 0; // Initialise
        signals[interactor_resite][ExperimentNo + offset] = 1;
    }
    else{ // if inserted before
        (signals[interactor_resite][ExperimentNo + offset]) = signals[interactor_resite][ExperimentNo + offset] + 1.0;
    }
//    cout << interactor_resite << "  " << signals[interactor_resite][ExperimentNo + offset] << endl;
}


