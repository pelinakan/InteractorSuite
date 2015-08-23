struct CaptureProbes{
	string chr;
    string name;
    string probeid;
    string target_id;
	int start;
	int end;
    int closestTSS_index;
	int annotated;// 0 = not annotated, 1 = annotated with a promoter, 2 = annotated with a cad_gwas, 3 = annotated with a negative control
	bool conflicting_annotations; // if a probe is associated with more than one feature
    
    FeatureStruct feature;
};

class ProbeSet{
    friend class PromoterClass;
public:
	vector <int> ChrRowStartIndexes;
	vector <int> ChrRowEndIndexes;
	vector <string> ChrNames;
	vector <CaptureProbes> Probes;
    
	void ReadProbeCoordinates();
	int AssociateReadwithProbes(string,int,int, PromoterClass);
    void AnnotateProbeswithPromoters(PromoterClass);
private:
	void GetProbeFeats(stringstream&, CaptureProbes&);
    int FindClosestTranscriptTSS(int, vector<int>);
};

void ProbeSet::GetProbeFeats(stringstream& line, CaptureProbes& t){
    
    string field;
    getline(line,t.chr,'\t');
    getline(line,field,'\t');
    t.start = atoi(field.c_str());
    getline(line,field,'\t');
    t.end = atoi(field.c_str());
    getline(line,t.target_id,'\t'); //target ID
    getline(line,t.probeid,'\t'); //probe ID
    getline(line,t.name,'\t'); //feature name
    
    if (t.target_id.find("refseq") != -1){
        t.annotated = 1;
        
    }
    else{
        if (t.target_id.find("gwas") != -1)
            t.annotated = 2;
        else
            t.annotated = 3;
    }
    
//    cout << t.target_id << "  " << t.annotated << endl;
}

void ProbeSet::ReadProbeCoordinates(){

	string filename;
	filename.append(dirname);
	filename.append("AgilentProbes_HG19_Design0415_all.txt");
	ifstream probefile(filename.c_str());

	long int index=0;
	string pline, chr1, chr2, temp;
    
    getline(probefile, temp);
    
    cout << temp;
	CaptureProbes tempprobe;

	getline(probefile,pline);
	stringstream probeline ( pline ); 
	GetProbeFeats(probeline,tempprobe);
	chr1=tempprobe.chr;
	ChrRowStartIndexes.push_back(index);
	ChrNames.push_back(chr1);
	Probes.push_back(tempprobe);
	++index;
	do{
		getline(probefile,pline);
        stringstream probeline ( pline );
		GetProbeFeats(probeline,tempprobe);
		chr2=tempprobe.chr;
		tempprobe.annotated = 0;
		tempprobe.conflicting_annotations = 0;
		Probes.push_back(tempprobe);
		++index;
		while(chr1 == chr2 && pline != ""){
			getline(probefile,pline);
            stringstream probeline ( pline );
			GetProbeFeats(probeline,tempprobe);
			chr2 = tempprobe.chr;
			Probes.push_back(tempprobe);
        //    cout << index << "   " << Probes[index].chr << "   " << Probes[index].start << "  " << Probes[index].target_id << "  " << Probes[index].annotated << "  " << endl;
			++index;
		}
     //   cout << "Probes on  " << chr1 << "  read" << endl;
		ChrRowEndIndexes.push_back(index-2);
		if(pline != ""){
			ChrRowStartIndexes.push_back(index-1);
			ChrNames.push_back(chr2);
		}
		chr1 = tempprobe.chr;
	//	cout << "Chr Index   " << index+1 << "   " << chr2 << endl;
	}while(!probefile.eof()&& pline!="" && pline != "END");
	Probes.pop_back();

	cout << Probes.size() <<  "   Probe Coordinates Read " << endl;
    
}


int ProbeSet::FindClosestTranscriptTSS(int probe_coord, vector<int> isoformcoords){
    
    int smallest_distance = 0, dist = 0, index = 0;
    smallest_distance = abs(probe_coord - isoformcoords[0]);
    for (int i = 0; i < isoformcoords.size(); ++i) {
        dist = abs(probe_coord - isoformcoords[i]);
        if (dist <= smallest_distance) {
            smallest_distance = dist;
            index = i;
        }
    }
    return index;
}

void ProbeSet::AnnotateProbeswithPromoters(PromoterClass proms){
    
    string pchr;
    vector<CaptureProbes>::iterator it;
    for(it = Probes.begin(); it != Probes.end(); ++it){
        if(it->annotated == 1){ //refseq
            pchr = it->chr;
            if(proms.GeneMap[pchr].Genes.find(it->name) != proms.GeneMap[pchr].Genes.end()){
                it->feature.chr = proms.GeneMap[pchr].Genes[it->name].chr;
                it->feature.RefSeqName = proms.GeneMap[pchr].Genes[it->name].RefSeqName;
                it->feature.strand = proms.GeneMap[pchr].Genes[it->name].strand;
                for (int j = 0; j < proms.GeneMap[pchr].Genes[it->name].isoformpromotercoords.size(); ++j) {
                    it->feature.TranscriptNames[j] = proms.GeneMap[pchr].Genes[it->name].TranscriptNames[j];
                    it->feature.isoformpromotercoords[j] = proms.GeneMap[pchr].Genes[it->name].isoformpromotercoords[j];
                }
                it->closestTSS_index = FindClosestTranscriptTSS(it->end, proms.GeneMap[pchr].Genes[it->name].isoformpromotercoords);
            }
            else
                it->annotated = 0;
        }
    }
}

int ProbeSet::AssociateReadwithProbes(string pchr, int readstart,int readend, PromoterClass proms){
int probe_index = -1, startsearch, endsearch;

startsearch = 0;
endsearch = -1;
    int i;
for(i = 0; i < ChrNames.size();++i){
	if(pchr == ChrNames[i]){
		startsearch = ChrRowStartIndexes[i];
		endsearch = ChrRowEndIndexes[i];
		break;
	}
}
 //   cout << startsearch << "  " << endsearch << endl;
for(probe_index = startsearch; probe_index <= endsearch; ++probe_index){
  //  cout << Probes[probe_index].start << "  " << readstart << "   " << Probes[probe_index].end << "  " << readend << "  " << probe_index << endl;
    if (CheckFragment_ifContainedwithinInterval((Probes[probe_index].start - padding),(Probes[probe_index].end + padding), readstart, readend)){
   //  cout << " on probe" << Probes[probe_index].start << "  " << readstart << "   " << Probes[probe_index].end << "  " << readend << "  " << probe_index << endl ;
        return probe_index;
    }
}
return -1;
}



