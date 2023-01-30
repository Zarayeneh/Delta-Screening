# README #
   
### How to get temporal communities? ###

* cd LifeSpan
    1. ***Input communities:*** community files from different time steps (community_t0, community_t1, ..., community_tn), each community file is a text file in which each line represents vertices within a commnity.
    2. ***Similarity threshold:***  first_threshold, as discribed in the paper, is ρ (core threshold)  for the percentage of common vertices of the current community and the core community, and the second_threshold is α (predecessor threshold) for the percentage of common vertices of the current community and the previous community.
    2. ***Run code:*** 
         - ```g++ -std=c++11 lifeCycle.cpp -o lifeCycle ``` 
         - ```./lifeCycle [-n num_timesteps] [-p first_threshold] [-q second_threshold] [-c name_input]" ``` 
    3. ***Output files:*** The code will output two files (LW.txt and out.txt) the first file contains the information of temporal communities (community id, start time, end time, average size, size of the community before it dies, and core community id) and the second file contains all information of the first file as well as the vertices of the core community. 

    #### Run the code with the provided examples ####
    Example to see the temporal communities is provided in folder LifeSpan/Example where each file is a community structure for a specific time step (30 time steps is chosen).
      
      - ```g++ -std=c++11 lifeCycle.cpp -o lifeCycle ``` 
      - ```./lifeCycle -n 30 -p 0.5 -q 0.5 -c comm ```
    
    