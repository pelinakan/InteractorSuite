/*
Background interaction frequency is calculated 
*/
class DetermineBackgroundLevels{
    friend class ProbeSet;
public:
	BG_signals bglevels;
	void CalculateMeanandStdRegress(ProbeSet&,string, int);
    void ReadBackgroundValues(string);

private:
	void LinearRegression(int,double*,double*,string);

};
void DetermineBackgroundLevels::CalculateMeanandStdRegress(ProbeSet& ngs, string BaseFileName, int ExperimentNo){


	boost::unordered::unordered_map< int, int > nofentries_perBin; // required for calculating mean
	boost::unordered::unordered_map< int, int > signal_square; // required for calculating stdev
	int distance;

	for(int i = 0; i < ngs.Probes.size(); i++){
        if (ngs.Probes[i].annotated == 3) {
            boost::unordered::unordered_map< int, int* >::const_iterator iter;
            for (iter = ngs.Probes[i].feature.Signals.signals.begin(); iter != ngs.Probes[i].feature.Signals.signals.end(); ++iter){
                distance = iter->first - ngs.Probes[i].end;
                int bin = abs(distance) / BinSize;
                if(bglevels.mean.find(bin) == bglevels.mean.end())
                    bglevels.mean[bin] = iter->second[ExperimentNo];
                else
                    bglevels.mean[bin] = bglevels.mean[bin] + iter->second[ExperimentNo];
                if(nofentries_perBin.find(bin) == nofentries_perBin.end()){
                    nofentries_perBin[bin] = 1;
                    signal_square[bin] = (iter->second[ExperimentNo])*(iter->second[ExperimentNo]);
                }
                else{
                    nofentries_perBin[bin] = nofentries_perBin[bin] + 1;
                    signal_square[bin] = (signal_square[bin] + ((iter->second[ExperimentNo])*(iter->second[ExperimentNo])));
                }
            }
        }
    }
// Calculate Mean and stdev
	boost::unordered::unordered_map< int, double>::iterator it; // iterator for bin signals
	for (it = bglevels.mean.begin(); it != bglevels.mean.end(); ++it){
		it->second = it->second / (double(nofentries_perBin[it->first])); // Mean for that bin
		double mean_square = ((it->second)*(it->second));
		double signalsquare_mean = double(signal_square[it->first] / double(nofentries_perBin[it->first]));
		bglevels.stdev[it->first] = sqrt((signalsquare_mean - mean_square));
	}

    string FileName;
	FileName.append(BaseFileName);
	FileName.append("BackgroundLevels.txt");
	ofstream outf(FileName.c_str());
	
    for (it = bglevels.mean.begin(); it != bglevels.mean.end(); ++it){
		outf << it->first << '\t' << it->second << '\t' << bglevels.stdev[it->first] << endl;
	}
}

void DetermineBackgroundLevels::ReadBackgroundValues(string BaseFileName){
    
    string FileName;
    FileName.append(BaseFileName);
    FileName.append("BackgroundLevels.txt");
    ifstream infile(FileName.c_str());
    
    int bin;
    double mean, stdev;
    while (!infile.eof()) {
        infile >> bin >> mean >> stdev;
        bglevels.mean[bin] = mean;
        bglevels.stdev[bin] = stdev;
    }
    cout << "Background Values Read" << endl;
}
void DetermineBackgroundLevels::LinearRegression(int n, double* x, double* y,string whichstream){

double intercept; //Power law y = (a) * (x^b)
int direction;

Maths::Regression::Linear A(n, x, y);


if (whichstream == "upstream")
	direction = 1;
else
	direction = 0;

	bglevels.b[direction] = A.getSlope();
    cout << "    Slope = " << bglevels.b[direction] << endl;
	
	intercept = A.getIntercept();
	cout << "Intercept = " << intercept << endl << endl;
   
	bglevels.a[direction] = pow(2,intercept);
	
	cout << "Regression coefficient = " << A.getCoefficient() << endl;

}


