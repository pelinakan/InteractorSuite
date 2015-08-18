//
//  SupplementaryFunctions.h
//  HiCapAnalysis
//
//  Created by Pelin Sahlen on 06/03/2015.
//  Copyright (c) 2015 Pelin Sahlen. All rights reserved.
//

#ifndef HiCapAnalysis_SupplementaryFunctions_h
#define HiCapAnalysis_SupplementaryFunctions_h

bool CheckPoint_ifContainedwithinInterval(int interval_start, int interval_end, int point){
    
    if (interval_end <= interval_start)
        cerr << "Given interval has a negative or zero length, aborting" << endl;
    
    if (point >= interval_start && point <= interval_end)
        return 1; // point contained within interval given, boundaries included within the interval
    else
        return 0;
    
}

bool CheckFragment_ifContainedwithinInterval(int interval_start, int interval_end, int frag_start, int frag_end){ //Interval coords first, then frag coords. Either frag_start or frag_end being in the interval is enough
    
    bool contained = 0;
    
    if (frag_end <= frag_start)
        cerr << "Given fragment has a negative or zero length, aborting" << endl;
    
   // cout << interval_start << "   " << interval_end << "   " << frag_start << "   " << frag_end << endl;
    
    if((CheckPoint_ifContainedwithinInterval(interval_start,interval_end,frag_start))){//If frag_start contained within interval
       contained = 1;
    }
    else{ // if frag_start not contained within interval, check if frag_end contained within interval
       if(CheckPoint_ifContainedwithinInterval(interval_start,interval_end,frag_end))
        contained = 1;
    }
    
    return contained;
}



#endif
