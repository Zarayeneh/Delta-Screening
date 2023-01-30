#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <typeinfo>
#include <sstream>
#include <unordered_set>
#include <algorithm>
#include <stdlib.h>     /* atoi */

using namespace std;

int N = 11;
float par1 = .5;
float par2 = .3;
char *filename = NULL;
unordered_map<int, int> cs;
class DisjSets
{
    unordered_map<int, int> parent;

public:

    // perform MakeSet operation
    void makeSet(vector<int> const &universe)
    {
        // create `n` disjoint sets (one for each item)
        for (int i: universe) {
            parent[i] = i;
        }
    }

    // Find the root of the set in which element `k` belongs
    int Find(int k)
    {
        // if `k` is root
        if (parent[k] == k) {
            return k;
        }

        // recur for the parent until we find the root
        return Find(parent[k]);
    }

    // Perform Union of two subsets
    void Union(int a, int b)
    {
        // find the root of the sets in which elements
        // `x` and `y` belongs
        int x = Find(a);
        int y = Find(b);

        parent[x] = y;
    }
};
//after making the dijoint data structure, now we want to see which nodes are in the same cluster and mark them
void make_cluster(DisjSets &dis,vector<int> v,vector<int> D)
{
    int r=0;
    for(int i=0;i<v.size();i++)
    {
        for(int j=r+1;j<r+v[i];j++)
        {
           dis.Union(D[r],D[j]);
        }
        r=r+v[i];
    }
 }

void print(vector<int> const &universe, DisjSets &dis) {
   for (int i : universe)
        cout << dis.Find(i) << " ";
   cout << '\n';
}

void print_vec(vector<int>& v)
{
    for(int i=0;i<v.size();i++)
        cout<<v[i]<< " ";
    cout<<endl;
}

vector<int> intersection(vector<int>& nums1, vector<int>& nums2) {
        unordered_set<int> m(nums1.begin(), nums1.end());
        vector<int> res;
        for (auto a : nums2)
            if (m.count(a)) {
                res.push_back(a);
                m.erase(a);
            }
        return res;
    }


void
usage(char *prog_name, const char *more) {
  cerr << more;
  cerr << "usage: " << prog_name << " [-n num_timesteps] [-p first_threshold] [-q second_threshold] [-c name_input]" << endl;
  cerr << "-n int\tshows number of time steps" << endl;
  cerr << "-p double\tshows the fisrt threshold, the percentage of common elements between current and previous time steps' communities" << endl;
  cerr << "-q double\tshows the second threshold, the percentage of common elements between current and core communities " << endl;
  cerr << "-c string\tname of input files" << endl;
  cerr << "-h\t show this usage message" << endl;

  exit(0);
}

void
parse_args(int argc, char **argv)
{
  if (argc<2)
    usage(argv[0], "Bad arguments number\n");

  for (int i = 1; i < argc; i++)
  {
        if(argv[i][0] == '-')
        {
          switch(argv[i][1])
          {
          case 'n':
              N = atoi(argv[i+1]);
              i++;
              break;
          case 'p':
              par1 = atof(argv[i+1]);
              i++;
              break;
          case 'q':
              par2 = atof(argv[i+1]);
              i++;
              break;
          case 'c':
              filename = argv[i+1];
              i++;
              break;
          case 'h':
              usage(argv[0], "");
              break;
          default:
            usage(argv[0], "Unknown option\n");
          }
        }
   }
}





int main(int argc, char** argv)
{
   //parse_args(argc, argv);;
   vector<int> comm_members0;
   vector<int> comm_members;
   vector<int> comm_size0;
   vector<int> comm_size;

   cout<<"k: "<< 0<<endl;
   //read the first input file containing nodes within each community
   string infile = "input0";
   cout<< "infile: "<<infile<<endl;
   ifstream myfile(infile);
   string line;
   while(getline(myfile,line))
   {
       stringstream ss(line);
       vector<int> v;
       int i;
       while(ss>>i)
       {
           v.push_back(i);
       }
       comm_size0.push_back(v.size());
       comm_members0.insert(comm_members0.end(), v.begin(), v.end());
   }
   
   myfile.close();
   DisjSets dis; //initialize DisjointSet class
   dis.makeSet(comm_members0);
   make_cluster(dis,comm_size0,comm_members0);//mark communities at time step t0
   //key: representative of core communities, value: first vector of t_start, t_end, AVG_size of previous life cycles, the size of the current comm and the second vector is members of the core community
   vector<unordered_map<int, pair<vector<int>,vector<int>> > > p(N);
   int r=0;
   for(int i=0;i<comm_size0.size();i++)//for the initial time step (t0) we consider all comm>=10 as core comm
   {
        vector<int> v(comm_members0.begin()+r, comm_members0.begin()+r+comm_size0[i]);
        vector<int> vp{0,0,comm_size0[i],comm_size0[i],dis.Find(comm_members0[r])};
        p[0].insert(make_pair(dis.Find(comm_members0[r]),make_pair(vp,v)));//
        r=r+comm_size0[i];
   }
   comm_size0.shrink_to_fit();
   comm_members0.shrink_to_fit();
   cout<<"L/W result: "<<endl;
   vector<int> vec_nodes;
   vec_nodes = comm_members0;
   /*for each time step find the community life cycle*/
   for(int k=1;k<N;k++)//N totall number of time steps
   {
       cout<<"k: "<<k<<endl;
       //read the input file containing nodes within each community commk
       string comm = "input"+to_string(k);
       ifstream myfile(comm);
       string line;
       while(getline(myfile,line))
       {
           stringstream ss(line);
           vector<int> v;
           int i;
           while(ss>>i)
            v.push_back(i);
            comm_size.push_back(v.size());
            comm_members.insert(comm_members.end(), v.begin(), v.end());
       }
       
       myfile.close();

       //for each community at the current time step this hash table shows key: rep comm of the prev comm, value: #common members
      // unordered_map<int,int> count;delete DisjSet;
       //key: keep the rep of the previous timestep satisfying the thresholds value: a vector of size two the first element is index of start of comm,index of the end of comm
       unordered_map<int,vector<int>> T;
       int r=0;
       for(int i=0;i<comm_size.size();i++)//for each comm in t1 we want to find the count of in common elements with prev timestep
       {
            unordered_map<int,int> count;
            for(int j=r;j<r+comm_size[i];j++)
            {
               std::vector<int>::iterator it = std::find(vec_nodes.begin(), vec_nodes.end(),comm_members[j]);//see if it is not a new node and exist at the p timestep
               if (it != vec_nodes.end())
               {

                   if(count.find(dis.Find(comm_members[j])) != count.end() )//if it is not the first time to see it increase the count by one
                   {
                       count[dis.Find(comm_members[j])]++;//
                   }
                   else
                   {
                       count.insert(make_pair(dis.Find(comm_members[j]),1));//if it is the first time to see insert one
                   }
               }
               else
               {
                   continue;//if it is a new node go to the next node
               }
            }

            for(auto c:count) 
            {//count has all common node counts with previous time step communities
                //if (intersection of current and previous comm)/(size of curr comm) && (intersection of current and previous comm)/(size of prev comm)  is mor than par1
                if(((double)c.second/(double)(p[k-1][c.first].first[3])>= par1 && ((double)c.second/(double)(comm_size[i])>= par1 )))//in common should be more than 70% we can consider jacard index but threshhold should be more than .54?
                {
                    vector<int> v(comm_members.begin()+r, comm_members.begin()+r+comm_size[i]);//members of the current comm
                    //if intersection with the root comm is mor than par2
                    int rep_inter = p[k-1][c.first].first[4];
                    int time_step = p[k-1][c.first].first[0];
                    vector<int> res = intersection(p[time_step][rep_inter].second,v);
                    if( (double)(res.size())/(double)(min((int)(p[time_step][rep_inter].first[3]),(int)comm_size[i])) >= par2)
                    {
                        //T: key:rep of previous comm, value: {index of start of comm,index of the end of comm}
                        vector<int> T_sec{c.first,comm_size[i],(p[k-1][c.first]).first[4]};
                        T[r]= T_sec;
                    }
                    else
                    {
                        vector<int> T_sec{-1,comm_size[i]};
                        T.insert(make_pair(r,T_sec));
                    }
                }
                else
                {
                    vector<int> T_sec{-1,comm_size[i]};
                    T.insert(make_pair(r,T_sec));
                }

            }
            r=r+comm_size[i];
       }

       dis.makeSet(comm_members);
       make_cluster(dis,comm_size,comm_members);//mark communities at time step t_1
       int rr=0;
       for(int ii=0;ii<comm_size.size();ii++)//for the initial time step (t0) we consider all comm>=10 as core comm
       {
               int start,ave_size,rep_root;
               vector<int> v1(comm_members.begin()+rr, comm_members.begin()+rr+comm_size[ii]);
               if(T.find(rr)!=T.end())
               {
                   if(T[rr][0] == -1)
                   {
                       start = k;
                       ave_size = comm_size[ii];
                       rep_root = dis.Find(comm_members[rr]);
                   }
                   else
                   {
                       vector<int> prev=p[k-1][T[rr][0]].first;
                       start = prev[0];
                       ave_size = (double)( (prev[2]* (prev[1]-prev[0]+1)) + T[rr][1] )/(double)(k-prev[0]+1);
                       rep_root = T[rr][2];
                   }
               }
               vector<int> vp1{start,k,ave_size,comm_size[ii],rep_root};//T[rr][1] is common members
               p[k].insert(make_pair(dis.Find(comm_members[rr]),make_pair(vp1,v1)));
               rr=rr+comm_size[ii];
       }


       vec_nodes = comm_members;
       for(auto m:p[k])
       {
            if((m.second.first[3]>=10) && (m.second.first[1] > m.second.first[0]))
                cout<<m.first<<"  " <<m.second.first[0] <<"  "<<m.second.first[1] <<"  " <<m.second.first[2]<<"  " <<m.second.first[3]<<"  " <<m.second.first[4]<<endl;
       }

       comm_size.clear();
       vector<int>().swap(comm_size);
       comm_members.clear();
       vector<int>().swap(comm_members);
       unordered_map<int,vector<int>>().swap(T);
   }
     //H/W result which shows the L/W along with the members of the core community
    ofstream out("out.txt", ios::app);
    for(int k=0;k<p.size();k++)
    {
       for(auto mp: p[k])
       {
           if(mp.second.first[0] != mp.second.first[1] && mp.second.first[3]>=10)// && mp.second.first[1] > mp.second.first[0])
           {
               out<<mp.first<<"\t" <<mp.second.first[0] <<"\t" <<mp.second.first[1]
               <<"\t" <<mp.second.first[2]<<"\t" <<mp.second.first[3]<<"\t";
               for(int j=0;j<p[mp.second.first[0]][mp.second.first[4]].second.size();j++)
                   out<<p[mp.second.first[0]][mp.second.first[4]].second[j]<<" ";
               out<<"\n";
           }
       }
    }
    out.close();
   return 0;
 }

