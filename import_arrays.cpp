


void import_basic(int* Default0,
                  int* Default,
                  int* Renew,
                  int* Policym,
                  int* Policya,
                  int* Policyh,
                  int* Policyr,
                  int* Policyl,
                  double* Puncond,
                  double* Pcond){


	// Default
	std::ifstream theFile ("matrices/Default1tr.txt");
  
    std::string line;
    std::vector<std::vector<std::string> > values;
    int it=0;
    std::string line_value;
    std::vector<std::string> line_values;
    std::stringstream ss;
    while(std::getline(theFile, line)){
      ss<<line;
      
      while(std::getline(ss, line_value, ',')){
        
        line_values.push_back(line_value);
        Default[it] = (int)::atof(line_value.c_str());
      }
      values.push_back(line_values);
      it=it+1;
      line_value.clear();
      line_values.clear();
      ss.clear();
    }


    std::ifstream theFile0 ("matrices/Default0tr.txt");
  
    it=0;
    while(std::getline(theFile0, line)){
      ss<<line;
      
      while(std::getline(ss, line_value, ',')){
        
        line_values.push_back(line_value);
        Default0[it] = (int)::atof(line_value.c_str());
      }
      values.push_back(line_values);
      it=it+1;
      line_value.clear();
      line_values.clear();
      ss.clear();
    }



    // Policym
    std::ifstream theFile2 ("matrices/Policym1tr.txt");
  
    it=0;
    while(std::getline(theFile2, line)){
      ss<<line;
      
      while(std::getline(ss, line_value, ',')){
        
        line_values.push_back(line_value);
        Policym[it] = (int)::atof(line_value.c_str());
      }
      values.push_back(line_values);
      it=it+1;
      line_value.clear();
      line_values.clear();
      ss.clear();
    }



    // Policya
    std::ifstream theFile3 ("matrices/Policya1tr.txt");
  
    it=0;
    while(std::getline(theFile3, line)){
      ss<<line;
      
      while(std::getline(ss, line_value, ',')){
        
        line_values.push_back(line_value);
        Policya[it] = (int)::atof(line_value.c_str());
      }
      values.push_back(line_values);
      it=it+1;
      line_value.clear();
      line_values.clear();
      ss.clear();
    }



    // Policyh
    std::ifstream theFile4 ("matrices/Policyh1tr.txt");
  
    it=0;
    while(std::getline(theFile4, line)){
      ss<<line;
      
      while(std::getline(ss, line_value, ',')){
        
        line_values.push_back(line_value);
        Policyh[it] = (int)::atof(line_value.c_str());
      }
      values.push_back(line_values);
      it=it+1;
      line_value.clear();
      line_values.clear();
      ss.clear();
    }



    // Policyl
    std::ifstream theFile5 ("matrices/Policyl1tr.txt");
  
    it=0;
    while(std::getline(theFile5, line)){
      ss<<line;
      
      while(std::getline(ss, line_value, ',')){
        
        line_values.push_back(line_value);
        Policyl[it] = (int)::atof(line_value.c_str());
      }
      values.push_back(line_values);
      it=it+1;
      line_value.clear();
      line_values.clear();
      ss.clear();
    }



    // Puncond
    std::ifstream theFile6 ("matrices/Puncond1tr.txt");
  
    it=0;
    while(std::getline(theFile6, line)){
      ss<<line;
      
      while(std::getline(ss, line_value, ',')){
        
        line_values.push_back(line_value);
        Puncond[it] = (double)::atof(line_value.c_str());
      }
      values.push_back(line_values);
      it=it+1;
      line_value.clear();
      line_values.clear();
      ss.clear();
    }




    // Pcond
    std::ifstream theFile7 ("matrices/Pcond1tr.txt");
  
    it=0;
    while(std::getline(theFile7, line)){
      ss<<line;
      
      while(std::getline(ss, line_value, ',')){
        
        line_values.push_back(line_value);
        Pcond[it] = (double)::atof(line_value.c_str());
      }
      values.push_back(line_values);
      it=it+1;
      line_value.clear();
      line_values.clear();
      ss.clear();
    }

}