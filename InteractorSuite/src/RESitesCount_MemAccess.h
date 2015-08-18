class RESitesClass{
	friend class PromoterClass;
public :
	int span;
	vector <string> chr_names;
	boost::unordered::unordered_map< string, int > chroffsets_indexfile;
	boost::unordered::unordered_map<string, int > chr_starts; // first RE site
	boost::unordered::unordered_map<string, int > chr_ends; // last RE site
	vector < REindexes > indexes;
	vector < int > posvector;
	void InitialiseVars();
	bool GettheREPositions(std::string, int, int*);
	int GetRESitesCount(string, int , int);
	void CleanClass(void);
private:
};
void RESitesClass::InitialiseVars(){

    string s;
	s.append(dirname);
	s.append("Digest_Human_GRCh37_MboI_None_13-48-29_04-08-2015.txt");
	ifstream RESitesf(s.c_str());
	s.clear();

	string chrname, chrp, temp;
	int pos, chrstart; 
	int startpos;
	int span = 1000000; // Window size
//For indexing
	cout << "Initialising RE site Class..." << endl;

bool f=0;

getline(RESitesf,temp); //get the header row1
getline(RESitesf,temp); //get the header row
RESitesf >> chrp >> temp >> pos >> temp >> temp >> temp >> temp; // Read the first line
chrname = chrp; // get the chr name outside the loop
chr_names.push_back(chrp);
chrstart = (pos ); // chromosome start
indexes.push_back(REindexes());
chroffsets_indexfile[chrp] = (indexes.size()-1);
chr_starts[chrp] = chrstart;
while(!(RESitesf.eof())){
	int binstart = (pos);
	int binend = binstart + span;
	startpos = posvector.size();
	while( chrname == chrp && (pos >= binstart && pos <= binend)){
		posvector.push_back(pos );
		RESitesf >> chrp >> temp >> pos >> temp >> temp >> temp >> temp;
		if(RESitesf.eof()){
			f = 1; //End of file
			break;
		}
	}
	if(indexes.back().binend.empty())
		indexes.back().binstart.push_back(chrstart);
	else
		indexes.back().binstart.push_back((indexes.back().binend.back())+1);
	indexes.back().binend.push_back((posvector.back() + 1));
	indexes.back().offset.push_back(startpos);
	indexes.back().count.push_back((posvector.size()-startpos));
	if((chrname != chrp) || f ){
		chr_ends[chrname] = indexes.back().binend.back();
		indexes.push_back(REindexes());
		chroffsets_indexfile[chrp] = (indexes.size()-1);
        chrname = chrp;
		chrstart = (pos);
		chr_names.push_back(chrp);
		chr_starts[chrp] = chrstart;
	}
}
	cout << "RE site class initialised " << endl;


}
int RESitesClass::GetRESitesCount(string chr, int st, int end){

int bitcount = 0;
int REcount = 0;
int rightchr;
int starttosearch = 0;
int HalfClusterDist = (coreprom_upstream + coreprom_downstream)/2;

boost::unordered::unordered_map< string, int >::iterator it = chroffsets_indexfile.find(chr);
rightchr = it->second; // Get the right index vector
for(int i = 0; i < indexes[rightchr].binstart.size();++i){ // Iterate over the bins
	if (indexes[rightchr].binstart[i] <= st && indexes[rightchr].binend[i] >= st){ // If start is within a bin
		starttosearch = indexes[rightchr].offset[i]; // mark the index in the posvector to start to search
		bitcount = indexes[rightchr].count[i]; // this many elements of the posvector is contained within that bin
		if (((st + HalfClusterDist) > indexes[rightchr].binend[i])) // if end is included in the next bin
			bitcount += indexes[rightchr].count[i+1]; // mark the number of elements in the next bin
			break;
	}
}
for (int i = starttosearch; i < starttosearch + bitcount ; ++i){ //This is where search in the posvector starts
	while ((posvector[i] >= st && posvector[i] <= end)){ // If that RE is within the fragment
		++i; 
		++REcount; // count RE sites
	}
	if(REcount > 0)
		break;
}
//cout << chr << ":" << st << "-" << end << "   " << REcount << endl;
	return REcount;

}

bool RESitesClass::GettheREPositions(std::string chr, int pos, int* renums){ // Returns closest RE sites to a position
	
    int HalfClusterDist = 5000; //in case the pos is at the end of a bin
    int rightchr;
    int starttosearch = 0;
    int bitcount = 0;
    int REposprev, REposat, REposnext;


    boost::unordered::unordered_map< string, int >::iterator it = chroffsets_indexfile.find(chr);
    
    if(it == chroffsets_indexfile.end()){ // Chromosome is not in the list
        //	cout << chr << "  cannot found.." << endl;
        return 0;
    }
    rightchr = it->second; // Get the right index vector
    boost::unordered::unordered_map< string, int >::iterator its = chr_starts.find(chr);
    boost::unordered::unordered_map< string, int >::iterator ite = chr_ends.find(chr);
    
    if ((pos ) <= its->second || (pos ) >= ite->second){
  //      cout << rightchr << "  " << its->second << "  " << ite->second  << "   " << pos << endl;
        return 0;
    }
    for(int i = 0; i < indexes[rightchr].binstart.size();++i){ // Iterate over the bins
        if (indexes[rightchr].binstart[i] <= (pos) && indexes[rightchr].binend[i] >= (pos)){ // If start-HalfClusterDist is within a bin
            if(i > 0){
                starttosearch = indexes[rightchr].offset[i-1]; // mark the index in the posvector to start to search
                bitcount = indexes[rightchr].count[i-1]; // this many elements of the posvector is contained within that bin
                bitcount += indexes[rightchr].count[i];
            }
            else{
                starttosearch = indexes[rightchr].offset[i]; // mark the index in the posvector to start to search
                bitcount = indexes[rightchr].count[i]; // this many elements of the posvector is contained within that bin
            }
            if (((pos + HalfClusterDist) > indexes[rightchr].binend[i]) && (i+1 < indexes[rightchr].binstart.size())) // if end is included in the next bin
                bitcount += indexes[rightchr].count[i+1]; // mark the number of elements in the next bin
            break;
        }
    }
    for (int i = (starttosearch + 1); i < starttosearch + bitcount; i++){
        REposprev = posvector[i - 1];
        REposat   = posvector[i];
        REposnext = posvector[i + 1];
        while (REposat <  pos){
            REposprev = REposat;
            ++i;
            REposat = posvector[i];
            REposnext = posvector[i + 1];
        }
        if (REposat == pos) {
            renums[0] = REposprev;
            renums[1] = REposnext;
            break;
        }
        else{
            renums[0] = REposprev;
            renums[1] = REposat;
            break;
        }
    }
    return 1;
}


void RESitesClass::CleanClass(void){

	posvector.clear();
	chr_names.clear();
	indexes.clear();
	chroffsets_indexfile.clear();
}