#include "Tree/TreeMap.h"

TreeMap::TreeMap( const string& stree ){
    //find the distance associated with the root node (if exists)
    string::const_iterator finish = stree.end();
    --finish;
    parentDistance = -1.0;
    while ( (*finish != ')') && (finish!=stree.begin()) ){
        --finish;
    }
    ++finish;
    while( isspace(*finish)&&(finish!=stree.end()) ){
        ++finish;
    }
    //read root distance
    if( (finish!=stree.end()) && (*finish==':') ){
        string::const_iterator intermediate = finish;
        ++intermediate;
        //strip space
        while(isspace(*intermediate)){
            ++intermediate;
        }
        string::const_iterator startNumber = intermediate;
        while( ( *intermediate != ';' ) &&
               ( !isspace(*intermediate ) ) &&
               ( intermediate != stree.end() ) ){
            if ( !isdigit(*intermediate) && (*intermediate!='.') ){
                cerr << "error in the distance of the root node : " << endl
                     << string(++finish, stree.end()) << endl;
            }
            ++intermediate;
        }
        string number(startNumber,intermediate);
        parentDistance = atof( number.c_str() );
    }
    treeMapConstruct( stree.begin(), finish );
}

TreeMap::TreeMap( string::const_iterator start,
                  string::const_iterator finish ) {
    treeMapConstruct( start, finish );
}

void TreeMap::treeMapConstruct( string::const_iterator start,
                           string::const_iterator finish ) {

    int c = count( start, finish, '(' );
    if ( c != count( start, finish, ')') ){
        cerr << "unbalanced label : " << endl
             << string(start, finish) << endl;
        exit(EXIT_FAILURE);
    }
    if (start == finish){
        cerr << "empty subtree" << endl;
        exit(EXIT_FAILURE);
    }
    while ( isspace(*start) ){
        ++start;
        if (start == finish){
            cerr << "empty subtree" << endl;
            exit(EXIT_FAILURE);
        }
    }
    --finish;
    while ( isspace(*finish) ){
        --finish;
        if (start == finish){
            cerr << "empty subtree" << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    //if not a leaf
    if ( c != 0 ){
        if ( ( *start != '(' ) || ( *finish != ')' ) ){
            cerr << "wrong subtree : " << endl
                 << string(start, ++finish) << endl;
            exit(EXIT_FAILURE);
        }
        string::const_iterator intermediate = start;
        ++intermediate;
        int depth;
        while ( intermediate != finish ){
            string::const_iterator startChild = intermediate;
            depth = 0;
            //go to the first free ',' or ':'
            //ie the first ',' or ':' not in a subtree
            //or the end of the string
            while ( !( ( (*intermediate == ',') || (*intermediate == ':') ) &&
                       ( depth == 0 ) )&&
                    ( intermediate != finish ) ){
                if ( *intermediate == '(' ){
                    ++depth;
                }
                if ( *intermediate == ')' ){
                    --depth;
                }
                ++intermediate;
                if ( depth < 0 ){
                    cerr << "error in the subtree : " << endl
                         << string(start, ++finish) << endl;
                    exit(EXIT_FAILURE);
                }                
            }
            assert(depth == 0);
            vector< pair<TreeMap,double> >::iterator iter =
                childMap.insert( childMap.end(), pair<TreeMap,double>
                    (TreeMap( startChild, intermediate ), -1.0) );
            //read a distance if provided
            if (*intermediate == ':'){
                ++intermediate;
                //strip space
                while(isspace(*intermediate)){
                    ++intermediate;
                }
                string::const_iterator startNumber = intermediate;
                while( ( *intermediate != ',' ) &&
                       ( !isspace(*intermediate ) ) &&
                       ( intermediate != finish ) ){
                    if ( !isdigit(*intermediate) && (*intermediate!='.') ){
                        cerr << "error in a distance of the subtree : " << endl
                             << string(start, ++finish) << endl;
                        exit(EXIT_FAILURE);
                    }
                    ++intermediate;
                }
                string number(startNumber,intermediate);
                (*iter).second = atof( number.c_str() );
                (*iter).first.setParentDistance((*iter).second);
                //strip space
                while( ( *intermediate != ',' ) &&
                       ( intermediate != finish ) ){
                    ++intermediate;
                }
            }
            //if necessary skip the ','
            if ( *intermediate == ',' ){
                ++intermediate;
                if(intermediate==finish){
                    cerr << "error in the subtree : " << endl
                         << string(start, ++finish) << endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
        if ( childMap.size() < 1 ){
            cerr << "error empty subtree in : " << endl
                 << string(start, ++finish) << endl;
            exit(EXIT_FAILURE);
        }
    }
    //if a leaf
    else{
        ++finish;
        label = string(start, finish);
    }
}


void TreeMap::printMap( ostream& outputStream ) const {
    //if it's a leaf
    if (childMap.size() == 0){
        outputStream << label;
    }
    else{
        outputStream << '(';
        for ( vector< pair<TreeMap,double> >::const_iterator iter = childMap.begin();
              iter != childMap.end(); ++iter ){
            if (iter != childMap.begin() ){
                outputStream << ',';
            }
            (*iter).first.printMap( outputStream );
            if ( (*iter).second != 0 ){
                outputStream << ':' << (*iter).second;
            }
        }
        outputStream << ')';
    }
}
