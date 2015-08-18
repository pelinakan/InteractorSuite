struct temppars{
    string chr;
    int start;
    int end;
    vector< double > expression;
    string strand;
    string name;
    string tr_id;
    string probe_id;
};

struct GeneList_perChr{
    boost::unordered::unordered_map< string, FeatureStruct > Genes;
};

class PromoterClass{ //Probe Clusters Associated with a Promoter
public:
    
    boost::unordered::unordered_map< string, GeneList_perChr > GeneMap;
    temppars *tp;
    void InitialiseData(void);
    void ReadPromoterAnnotation(RESitesClass&, MappabilityClass&);
    bool AnnotatewithPromoters(string, int, int);
    int NofGenes;
    int NofTranscripts;
    
private:
    void GetTrFeats(stringstream&,temppars&);
    int FindLeftMostCoord(vector<int>);
    int FindRightMostCoord(vector<int>);
};

void PromoterClass::InitialiseData(void){
    
    tp = new temppars [2];
    
    cout << "Promoter Class Initialised" << endl;
    
}
void PromoterClass::GetTrFeats(stringstream &trx, temppars &tpars){
    
    string field;
    getline(trx,tpars.name,'\t');
    
    getline(trx,tpars.tr_id,'\t'); // tr id
    getline(trx,tpars.chr,'\t');
    getline(trx,tpars.strand,'\t');
    getline(trx,field,'\t');
    if(tpars.strand=="+"){
        tpars.start=atoi(field.c_str());
        getline(trx,field,'\t');
        tpars.end=atoi(field.c_str());
    }
    else{
        tpars.end=atoi(field.c_str());
        getline(trx,field,'\t');
        tpars.start=atoi(field.c_str());
    }
    //getline(trx,tpars.probe_id,'\t');
    getline(trx,field,'\t');
    getline(trx,field,'\t');
    getline(trx,field,'\t');
    
    getline(trx,field,'\t');
    tpars.expression.push_back(atof(field.c_str())); //mES
    getline(trx,field,'\t');
    tpars.expression.push_back(atof(field.c_str())); //XEN
    getline(trx,field,'\t');
    tpars.expression.push_back(atof(field.c_str())); //TS

}

int PromoterClass::FindLeftMostCoord(vector<int> coords){
    int leftmosttss;
    
    leftmosttss = coords[0];
    for (int i = 0; i < coords.size();++i){
        if (coords[i] <= leftmosttss )
            leftmosttss = coords[i];
    }
    return leftmosttss;
}

int PromoterClass::FindRightMostCoord(vector<int> coords){
    int rightmosttss;
    
    rightmosttss = coords[0];
    for (int i = 0; i < coords.size();++i){
        if (coords[i] >= rightmosttss )
            rightmosttss = coords[i];
    }
    return rightmosttss;
    
}

void PromoterClass::ReadPromoterAnnotation(RESitesClass& dpnIIsites, MappabilityClass& mapp)
{
    string temp,tr1, tr2;
    vector< int > isoformpromotercoords;
    vector< string > tr_ids;
    
    string RefSeqfilename;
    RefSeqfilename.append(dirname);
    RefSeqfilename.append("hg19.refseq.txt");
    ifstream RefSeq_file(RefSeqfilename.c_str());
    
    getline(RefSeq_file,temp); //get the header row

    getline(RefSeq_file,tr1);
    stringstream trx1 ( tr1 );
    GetTrFeats(trx1,tp[0]);
    isoformpromotercoords.push_back(tp[0].start);
    tr_ids.push_back(tp[0].tr_id);
    do{
        getline(RefSeq_file,tr2);
        stringstream trx2 ( tr2);
        GetTrFeats(trx2,tp[1]);
        while(tp[0].name == tp[1].name){
            ++NofTranscripts;
            isoformpromotercoords.push_back(tp[1].start);
            tr_ids.push_back(tp[1].tr_id);
            getline(RefSeq_file,tr2);
            stringstream trx2 ( tr2);
            GetTrFeats(trx2,tp[1]);
        };
        
        GeneMap[tp[0].chr].Genes[tp[0].name].chr = tp[0].chr;
        GeneMap[tp[0].chr].Genes[tp[0].name].strand = tp[0].strand;
        GeneMap[tp[0].chr].Genes[tp[0].name].RefSeqName = tp[0].name;
        for (int i = 0; i < isoformpromotercoords.size(); ++i) {
            GeneMap[tp[0].chr].Genes[tp[0].name].TranscriptNames.push_back(tr_ids[i]);
            GeneMap[tp[0].chr].Genes[tp[0].name].isoformpromotercoords.push_back(isoformpromotercoords[i]);
        }        
        isoformpromotercoords.clear();
        tr_ids.clear();
        //SWAP TP[0] and TP[1]
        tp[0].name.clear();
        tp[0].name.append(tp[1].name);
        tp[0].tr_id.clear();
        tp[0].tr_id.append(tp[1].tr_id);
        tp[0].chr.clear();
        tp[0].chr.append(tp[1].chr);
        tp[0].end = tp[1].end;
        tp[0].start = tp[1].start;
        tp[0].strand = tp[1].strand;
        
        ++NofGenes;
        
        if(NofGenes % 5000 == 0){
            cout << NofGenes << "    Genes Read" << endl;
            //break;
        }
    }while(tp[1].name!="END");
    
    cout << NofGenes << "  genes corresponding to " << NofTranscripts << "  transcripts are read from " << RefSeqfilename << endl;
}
bool PromoterClass::AnnotatewithPromoters(string chr, int pos_firstinpair, int len_firstinpair){
    
    bool pann = 0;
    boost::unordered::unordered_map< string, FeatureStruct >::iterator geneit;
    
    if(GeneMap.find(chr) != GeneMap.end()){
        for (geneit = GeneMap[chr].Genes.begin(); geneit != GeneMap[chr].Genes.end();++geneit){
            vector < int >::iterator coordit;
            for (coordit = geneit->second.isoformpromotercoords.begin(); coordit != geneit->second.isoformpromotercoords.end();++coordit) {
                if (CheckFragment_ifContainedwithinInterval((*coordit - coreprom_upstream), (*coordit + coreprom_downstream), pos_firstinpair, (pos_firstinpair + len_firstinpair))){ // If any base of the first read is contained within the core promoter
                    pann = 1; // Read is annotated with a promoter
                }
            }
        }
    }
    cout << " Annotated with prom " << pann << endl;
    return pann;
}






