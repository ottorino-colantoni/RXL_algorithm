/*
 * Auxiliary.cpp
 *
 *  Created on: 01 feb 2016
 *      Author: Mattia D'Emidio
 */

#include "Auxiliary.h"
#include "networkit/components/ConnectedComponents.hpp"
#include "networkit/components/StronglyConnectedComponents.hpp"
#include "networkit/generators/ErdosRenyiGenerator.hpp"
#include "networkit/generators/BarabasiAlbertGenerator.hpp"
#include <unordered_set>
#include <tuple>
#include <boost/algorithm/string.hpp>
#include <networkit/centrality/EstimateBetweenness.hpp>

Auxiliary::Auxiliary() {}

std::pair<std::vector<custom_node>,std::vector<custom_node>> Auxiliary::compute_ordering(NetworKit::Graph* graph,int ordering) {


	std::vector<custom_node> order;
	std::vector<custom_node> reverse_order;
	INFO("Computing Ordering");


	order.resize(graph->upperNodeIdBound(),NULL_NODE);
	reverse_order.resize(graph->upperNodeIdBound(),NULL_NODE);

	std::vector<std::pair<std::pair<custom_node, custom_node>, custom_node> > ranked_nodes(graph->upperNodeIdBound());

	std::cout<<std::flush;
	assert(ordering==0 || ordering==1 || ordering==2);

	if(ordering==2){
		std::cout<<"Manual order setting\n";
		int n = 0;
		for(size_t t=0;t<order.size();t++){
			std::cout<<"Insert vertex in position "<<t<<"\n";

			std::cin>>n;
			while(std::cin.fail()){
				std::cout<<"Insert (integer) vertex in position "<<t<<"\n";
			}
			order[t]=n;
			reverse_order[n]=t;

		}
		return std::make_pair(order,reverse_order);
	}
	if(ordering==0){

		INFO("BY DEGREE");
		graph->parallelForNodes([&] (custom_node i){
			assert(graph->hasNode(i));
			ranked_nodes[i] = graph->isDirected() ?
					std::make_pair(std::make_pair(graph->degreeIn(i)+graph->degreeOut(i), rand()), i) :
					std::make_pair(std::make_pair(graph->degree(i), rand()), i);

		});
	}
	else{

		INFO("BY APX BTW 2");
		const NetworKit::Graph& hand = *graph;
		std::cout<<round(std::pow(graph->upperNodeIdBound(),50.0/100.0)) <<" SAMPLES\n";
		NetworKit::EstimateBetweenness* rank = new NetworKit::EstimateBetweenness
				(hand,round(std::pow(graph->upperNodeIdBound(),50.0/100.0)) ,false,true);
		rank->run();

		graph->parallelForNodes([&] (custom_node i){
			if(std::isnan(rank->score(i))==1) //NaN removal
				ranked_nodes[i] = std::make_pair(std::make_pair(0, rand()), i);
			else
				ranked_nodes[i] = std::make_pair(std::make_pair(rank->score(i), rand()), i);

		});

		delete rank;


	}

	sort(ranked_nodes.begin(), ranked_nodes.end());
	reverse(ranked_nodes.begin(), ranked_nodes.end());

	for(size_t count = 0; count < ranked_nodes.size();count++){
		if(count>=graph->numberOfNodes())
			break;
		order[count] = ranked_nodes[count].second;
		reverse_order[order[count]] = count;
		if(count<10)
			std::cout<<"Position "<<count<<" of vertex "<<order[count]<<" rank "<<ranked_nodes[count].first.first<<std::endl;
	}
	return std::make_pair(order,reverse_order);
}

void Auxiliary::test_match(NetworKit::Graph* g,Labeling* l,int soglia){

	const NetworKit::Graph& G = *g;

	for (int var = 0; var < soglia; ++var) {
		custom_node i = G.randomNode();
		custom_node j = G.randomNode();
		if(!g->isWeighted()){
			NetworKit::BFS* dij = new NetworKit::BFS(G,i,true,false,j);
			dij->run();
			if(dij->distance(j)!=l->query(i,j)){
				std::cout<<dij->distance(j)<<" "<<l->query(i,j)<<"\n";
				throw new std::runtime_error("difference");
			}
			delete dij;
		}
		else{
			NetworKit::Dijkstra* dij = new NetworKit::Dijkstra(G,i,true,false,j);
			dij->run();
			if(dij->distance(j)!=l->query(i,j)){
				std::cout<<dij->distance(j)<<" "<<l->query(i,j)<<"\n";
				throw new std::runtime_error("difference");
			}
			delete dij;
		}
	}

}

void Auxiliary::store_mods(std::string fname,std::vector<NetworKit::GraphEvent>& mods){

	std::ofstream ofs(fname);

	if (!ofs)
		throw std::runtime_error("Error opening File ");

	int increm_ops = 0;
	int decrem_ops = 0;
	for(auto&el:mods){
		if(el.type==NetworKit::GraphEvent::EDGE_REMOVAL || (el.type==NetworKit::GraphEvent::EDGE_WEIGHT_UPDATE && el.w>0))
			decrem_ops++;
		else
			increm_ops++;
	}
	ofs<<boost::lexical_cast<std::string>(increm_ops)<<" "<<boost::lexical_cast<std::string>(decrem_ops)<<"\n";


	ProgressStream writer(mods.size());

	writer.label() << "Writing "<< mods.size()<<" mods to file ";
	for(auto&el:mods){
		switch (el.type) {
		case NetworKit::GraphEvent::EDGE_REMOVAL:{
			ofs<<"D "<<el.u<<" "<<el.v<<" "<<el.w<<"\n";
			break;
		}

		case NetworKit::GraphEvent::EDGE_ADDITION:{
			ofs<<"I "<<el.u<<" "<<el.v<<" "<<el.w<<"\n";
			break;
		}
		case NetworKit::GraphEvent::EDGE_WEIGHT_UPDATE:{
			if(el.w>0){
				ofs<<"WI "<<el.u<<" "<<el.v<<" "<<el.w<<"\n";
			}
			else{
				ofs<<"WD "<<el.u<<" "<<el.v<<" "<<-el.w<<"\n";
			}
			break;
		}

		default:
			break;
		}
		++writer;


	}


	ofs.close();

}

void Auxiliary::init_mods(std::vector<NetworKit::GraphEvent>& output_vector, int modT, int num, NetworKit::Graph* graph, double probability_of_weight_changes){


	/* initialize random seed: */
	MersenneTwister* randomGenerator;
	randomGenerator = new MersenneTwister();
	unsigned int actual_mod_type = modT;
	custom_node n1,n2;

	if(!graph->isWeighted())
		probability_of_weight_changes = 0;
	while(output_vector.size()<num){

		if(modT==2)//fully, flip coin
			actual_mod_type = randomGenerator->getRandomInteger()%2;

		NetworKit::GraphEvent event;

		if(actual_mod_type==0){ //decremental
			std::random_device rd;
		    std::mt19937 e2(rd());
		    std::uniform_real_distribution<> dist(0, 1);
			std::pair<custom_node,custom_node> edge = graph->randomEdge();
			n1 = edge.first;
			n2 = edge.second;
			event.u = n1;
			event.v = n2;
		    if(dist(e2)<=probability_of_weight_changes){
		    	event.type = NetworKit::GraphEvent::EDGE_WEIGHT_UPDATE;
				double actual_range = (double)(randomGenerator->getRandomInteger()%100 + 1)/100;
				assert(actual_range<=1);
				event.w = graph->weight(n1,n2)+round(actual_range*graph->weight(n1,n2));
				assert(event.w>graph->weight(n1,n2));
		    }
		    else{
		    	event.type = NetworKit::GraphEvent::EDGE_REMOVAL;
				event.w = -1;
		    }

		}
		else{
			//incremental
			std::random_device rd;
			std::mt19937 e2(rd());
			std::uniform_real_distribution<> dist(0, 1);
			if(dist(e2)<=probability_of_weight_changes){
				while(true){
					std::pair<custom_node,custom_node> edge = graph->randomEdge();
					n1 = edge.first;
					n2 = edge.second;
					event.u = n1;
					event.v = n2;
					event.type = NetworKit::GraphEvent::EDGE_WEIGHT_UPDATE;
					double actual_range = (double)(randomGenerator->getRandomInteger()%100 + 1)/100;
					event.w = graph->weight(n1,n2) - round(actual_range*graph->weight(n1,n2)+0.5);
					assert(event.w<graph->weight(n1,n2));

					if(round(actual_range*graph->weight(n1,n2))>0 && event.w>0){
						event.w = -event.w;
						break;
					}

				}
			}
			else{
				while(true){
					n1 = graph->randomNode();
					n2 = graph->randomNode();

					if(!graph->hasEdge(n1,n2) && n1!=n2)
						break;
				}
				event.u = n1;
				event.v = n2;
				event.type = NetworKit::GraphEvent::EDGE_ADDITION;
				std::pair<custom_node,custom_node> edge = graph->randomEdge();
				event.w = graph->weight(edge.first,edge.second);
			}

		}
		bool exists = false;
		for(auto& element : output_vector){
			if(element.u == event.u && element.v == event.v){
				exists = true;
				break;
			}
			if(!graph->isDirected()){
				if (element.u == event.v && element.v == event.u){
					exists = true;
					break;
				}
			}
		}
		if(exists)
			continue;
		else{
			output_vector.push_back(event);
		}
	}
	delete randomGenerator;
}

void Auxiliary::read_mods(std::string source, std::vector<NetworKit::GraphEvent>& output_vector) {

	std::ifstream ifs(source);
	if (!ifs)
		throw std::runtime_error("Error opening File ");

	std::cout<<"Reading mods from " << source << std::endl;


	int num_incremental = -1, num_decremental = -1;
	ifs >> num_incremental >> num_decremental;
	ProgressStream reader(num_incremental+num_decremental);
	reader.label() << "Reading "<< num_incremental+num_decremental << " Dynamic Operations ";
	custom_node n1,n2;
	std::string weight,type;
	while(true){
		ifs >> type >> n1 >> n2 >> weight;

		if(ifs.eof())
			break;

		++reader;
		assert(boost::lexical_cast<custom_weight>(weight)!=0);


		if(boost::iequals(type, "D"))
			output_vector.push_back(NetworKit::GraphEvent(NetworKit::GraphEvent::EDGE_REMOVAL,n1,n2,-1));
		else if(boost::iequals(type, "I")){
			//add with weight read from file
			output_vector.push_back(NetworKit::GraphEvent(NetworKit::GraphEvent::EDGE_ADDITION,n1,n2,boost::lexical_cast<custom_weight>(weight)));
		}
		else if(boost::iequals(type, "WD")){
			//weight decrease with negative values
			output_vector.push_back(NetworKit::GraphEvent(NetworKit::GraphEvent::EDGE_WEIGHT_UPDATE,n1,n2,-boost::lexical_cast<custom_weight>(weight)));
		}
		else if(boost::iequals(type, "WI")){
			output_vector.push_back(NetworKit::GraphEvent(NetworKit::GraphEvent::EDGE_WEIGHT_UPDATE,n1,n2,boost::lexical_cast<custom_weight>(weight)));
		}
		else
			throw new std::runtime_error("undefined weight");




	}
	ifs.close();



	int count_inc=0,count_decr=0;
	for(auto&el:output_vector)
		if(el.type==NetworKit::GraphEvent::EDGE_REMOVAL || (el.type==NetworKit::GraphEvent::EDGE_WEIGHT_UPDATE && el.w>0))
			count_decr++;
		else
			count_inc++;

	if(output_vector.size()!=num_incremental+num_decremental || count_decr!=num_decremental || count_inc!=num_incremental)
		throw std::runtime_error("Wrong mods reading");




}


bool Auxiliary::isInteger(const std::string & s){
   if(s.empty() || ((!isdigit(s[0])) && (s[0] != '-') && (s[0] != '+'))) return false ;

   char * p ;
   strtol(s.c_str(), &p, 10) ;

   return (*p == 0) ;
}
void Auxiliary::get_all(const boost::filesystem::path& root, const std::string& ext, std::vector<boost::filesystem::path>& ret){

    if(!boost::filesystem::exists(root) || !boost::filesystem::is_directory(root)) return;

    boost::filesystem::recursive_directory_iterator it(root);
    boost::filesystem::recursive_directory_iterator endit;

    while(it != endit){
        if(boost::filesystem::is_regular_file(*it) && it->path().extension() == ext) ret.push_back(it->path().filename());
        ++it;
    }
}
double Auxiliary::getQuartile(int type, std::vector<double>& v){

    std::sort(v.begin(),v.end());

    size_t sizeV = v.size();

    //determine quartiles
    size_t midSpeedUp = sizeV / 2;
    size_t qOne = midSpeedUp / 2;
    size_t qThree = (sizeV + midSpeedUp) / 2;
    switch(type){
      case 0:{
	return v[0];

      }
      case 1:{
	return v[qOne];
      }
      case 2:{
	return sizeV % 2 == 0 ? (v[midSpeedUp] + v[midSpeedUp-1]) /2 : v[midSpeedUp];
	}
      case 3:{
	return v[qThree];
	}
      case 4:{
	return v[v.size()-1];
	}
  default: throw std::runtime_error("unexpected behavior");


    }


}
void Auxiliary::box_plot_column(std::string data_folder,int pos) {
	std::vector<boost::filesystem::path> fileList;
	get_all(data_folder,".res",fileList);


	if(fileList.size()==0)
	  throw std::runtime_error("empty input folder");

	std::string outFolder = data_folder+"/output/";
	boost::filesystem::path dir(outFolder.c_str());
	if(!boost::filesystem::exists(dir))
	  boost::filesystem::create_directory(dir);

	std::string pFile = outFolder+"plotter.p";
	std::fstream ofs_plot;
	ofs_plot.open(pFile.c_str(),std::ios::out);
	ofs_plot<<"set terminal epslatex\n";
	ofs_plot<<"set out 'graph.tex'\n";
	ofs_plot<<"unset key\n";
	ofs_plot<<"set boxwidth 0.5 absolute\n";
	std::vector<std::string> labels;
	for(size_t it=0;it<fileList.size();it++){
	  boost::filesystem::path path_element = fileList[it];
	  std::cout<<"Insert label for file: "<<path_element.filename().string()<<" : ";
	  std::string l;
	  std::cin>>l,
	  labels.push_back(l);
	}
	ofs_plot<<"set xrange [-1:"<<2*labels.size()-1<<"]\n";
	ofs_plot<<"set logscale y\n";
	ofs_plot<<"set style line 1 lc rgb '#0060ad' lt 1 lw 1 pt 7 ps 1\n";
	ofs_plot<<"set style fill empty\n";
	for(size_t it=0;it<fileList.size();it++){
		  boost::filesystem::path path_element = fileList[it];
	      std::ifstream infile;
	      infile.open(path_element.string());
	      if(!infile)
	    	  throw std::runtime_error("error opening file");

	      std::vector<double> data;

	      std::string commands = " wc -l " + path_element.string() +" | awk '{print $1}'";
	      int n_lines = atoi(exec(commands.c_str()).c_str());
	      std::string line;
			INFO("WARNING CHECK POS<LENGTH OF ROW");

		  ProgressStream reader(n_lines-2);
		  reader.label() << "Processing file " << path_element.filename() << " with "<< n_lines-2 << " rows";
	      std::string temp;
	      getline(infile,line);
	      int counter = 0;
	      while (true){

				getline(infile,line);
				counter++;
				++reader;
				if(infile.eof() || counter == (n_lines-2)) //skip last two rows
				  break;



				std::stringstream ss(line);

				int col = 0;
				double local_data;
				while(col<=pos){
					ss >> temp;
					col++;
				}

				ss >> local_data;
				data.push_back(local_data);
		  }
//std::sort(data.begin(),data.end());
//for(auto& el: data){
//	std::cout<<el<<"\n";
//}

		  infile.close();
		  std::string filespeedUp = outFolder+path_element.filename().string()+".speedup";
		  std::fstream ofs;

		  ofs.open(filespeedUp.c_str(),std::ios::out);
		  ofs<<it*2<<" ";
		  for(int i=0;i<=4;i++)
			  ofs<<getQuartile(i,data)<<" ";
		  ofs<<"0.5 "<<labels[it]<<std::endl;
		  ofs.close();

		  if(it==0){
			  ofs_plot<<"plot '"<<path_element.filename().string()<<".speedup'"<<
					  " using 1:3:2:6:5:7:xticlabels(8) with candlesticks notitle whiskerbars ls 1, \\\n"
					  <<"'' using 1:4:4:4:4:7 with candlesticks lt -1 notitle, \\\n";
		  }

		  else if(it<fileList.size()-1){
			  ofs_plot<<" '"<<path_element.filename().string()<<".speedup'"<<
					  " using 1:3:2:6:5:7:xticlabels(8) with candlesticks notitle whiskerbars ls 1, \\\n"
					  <<"'' using 1:4:4:4:4:7 with candlesticks lt -1 notitle, \\\n";
		  }

		  else{
			  ofs_plot<<" '"<<path_element.filename().string()<<".speedup'"<<
					  " using 1:3:2:6:5:7:xticlabels(8) with candlesticks notitle whiskerbars ls 1, \\\n"
					  <<"'' using 1:4:4:4:4:7 with candlesticks lt -1 notitle";
		  }

		}

		ofs_plot.close();


}
void Auxiliary::process_data(std::string data_folder) {

	std::vector<boost::filesystem::path> fileList;
	get_all(data_folder,".res",fileList);


	if(fileList.size()==0)
	  throw std::runtime_error("empty input folder");

	std::string outputFolder = data_folder+"/output/";
	boost::filesystem::path dir(outputFolder.c_str());

	std::string plotterFile = outputFolder+"plotter.p";
	std::string statsVarieFile = outputFolder+"miscstats.dat";
	std::fstream ofs_plot,ofs_stats;
	ofs_plot.open(plotterFile.c_str(),std::ios::out);
	ofs_stats.open(statsVarieFile.c_str(),std::ios::out);
	ofs_plot<<"set terminal epslatex\n";
	ofs_plot<<"set out 'graph.tex'\n";
	ofs_plot<<"unset key\n";
	ofs_plot<<"set boxwidth 0.5 absolute\n";


	if(!boost::filesystem::exists(dir))
	  boost::filesystem::create_directory(dir);

	std::vector<std::string> names;
	for(size_t it=0;it<fileList.size();it++){
	  boost::filesystem::path path_element = fileList[it];
	  //std::size_t pos = path_element.filename().string().find(".txt");
	  std::size_t pos = 10;
	  names.push_back(path_element.filename().string().substr(0,pos));
	}
	ofs_plot<<"set xrange [-1:"<<2*names.size()-1<<"]\n";
	ofs_plot<<"set logscale y\n";
	ofs_plot<<"set style line 1 lc rgb '#0060ad' lt 1 lw 1 pt 7 ps 1\n";
	ofs_plot<<"set style fill empty\n";

	for(size_t it=0;it<fileList.size();it++){


		  boost::filesystem::path path_element = fileList[it];


	      std::ifstream infile;
	      infile.open(path_element.string());
	      if(!infile)
	    	  throw std::runtime_error("error opening file");

	      std::vector<double> timeFS,timeDyn,ratioTime;
	      std::vector<double> spaceFS,spaceDyn,ratioSpace;
	      std::vector<double> queryFS,queryDyn,ratioQuery;
	      std::vector<double> affected;
	      std::vector<double> hatA;

	      std::string commands = " wc -l " + path_element.string() +" | awk '{print $1}'";
	      int n_lines = atoi(exec(commands.c_str()).c_str());

	      std::string line;
	      ProgressStream reader(n_lines-2);
	      reader.label() << "Processing file " << path_element.filename() << " with "<< n_lines-2 << " rows";


	      double t1,t2,t3;
	      double s1,s2,s3;
	      double q1,q2,q3;
	      double aff,aff_x,aff_y;
	      std::string temp;

	      getline(infile,line);
	      int counter = 0;

	      while (true){

				getline(infile,line);
				counter++;
				++reader;
				if(infile.eof() || counter == (n_lines-2)) //skip last two rows
				  break;



				std::stringstream ss(line);
				ss >> temp >> t1 >> t2 >> t3 >> q1 >> q2 >> q3 >> s1 >> s2 >> s3 >> aff >> aff_x >> aff_y;

				//std::cout << temp << " "<< t1 << " " << t2 << " "<< t3 << " " <<aff_x << " "<<aff_y<<std::endl;
				timeFS.push_back(t1);

				timeDyn.push_back(t2);
				ratioTime.push_back(t3);
				spaceFS.push_back(s1);
				spaceDyn.push_back(s2);
				ratioSpace.push_back(s3);
				queryFS.push_back(q1);
				queryDyn.push_back(q2);
				ratioQuery.push_back(q3);
				affected.push_back(aff);
				if(aff_x>aff_y)
					hatA.push_back(aff_x);
				else
					hatA.push_back(aff_y);
	      }


	      infile.close();

	      std::string filespeedUp = outputFolder+path_element.filename().string()+".speedup";
	      std::fstream ofs;

	      ofs.open(filespeedUp.c_str(),std::ios::out);
	      ofs<<it*2<<" ";
	      for(int i=0;i<=4;i++)
	    	  ofs<<getQuartile(i,ratioTime)<<" ";
	      ofs<<"0.5 "<<names[it]<<std::endl;
	      ofs.close();

	      if(it==0){
	    	  ofs_stats<<"graph\t\t\tavgTime\t\tstdDev\t\tavgAff\t\tavgHatA"<<std::endl;
	    	  ofs_stats<<path_element.filename().string()<<"\t\t"<<avg(timeDyn)<<"\t\t"<<stddev(timeDyn)<<"\t\t"<<avg(affected)<<"\t\t"<<avg(hatA)<<std::endl;
	    	  ofs_plot<<"plot '"<<path_element.filename().string()<<".speedup'"<<
	    			  " using 1:3:2:6:5:7:xticlabels(8) with candlesticks notitle whiskerbars ls 1, \\\n"
	    			  <<"'' using 1:4:4:4:4:7 with candlesticks lt -1 notitle, \\\n";
	      }

	      else if(it<fileList.size()-1){
	    	  ofs_stats<<path_element.filename().string()<<"\t\t"<<avg(timeDyn)<<"\t\t"<<stddev(timeDyn)<<"\t\t"<<avg(affected)<<"\t\t"<<avg(hatA)<<std::endl;
	    	  ofs_plot<<" '"<<path_element.filename().string()<<".speedup'"<<
	    			  " using 1:3:2:6:5:7:xticlabels(8) with candlesticks notitle whiskerbars ls 1, \\\n"
	    			  <<"'' using 1:4:4:4:4:7 with candlesticks lt -1 notitle, \\\n";
	      }

	      else{
	    	  ofs_stats<<path_element.filename().string()<<"\t\t"<<avg(timeDyn)<<"\t\t"<<stddev(timeDyn)<<"\t\t"<<avg(affected)<<"\t\t"<<avg(hatA)<<std::endl;
	    	  ofs_plot<<" '"<<path_element.filename().string()<<".speedup'"<<
	    			  " using 1:3:2:6:5:7:xticlabels(8) with candlesticks notitle whiskerbars ls 1, \\\n"
	    			  <<"'' using 1:4:4:4:4:7 with candlesticks lt -1 notitle";
	      }

	    }
	    ofs_stats.close();
	    ofs_plot.close();


}

void Auxiliary::SCC(NetworKit::Graph* source_graph){

	std::cout<<"Original graph has "<<source_graph->numberOfNodes()<<" vertices and "<<source_graph->numberOfEdges()<<" edges"<<std::endl;

	NetworKit::Partition p;

	if(source_graph->isDirected()){
		NetworKit::StronglyConnectedComponents* scc = new NetworKit::StronglyConnectedComponents(*source_graph);
		scc->run();
		p = scc->getPartition();
		delete scc;
	}
	else{
		NetworKit::ConnectedComponents* conn = new NetworKit::ConnectedComponents(*source_graph);
		conn->run();
		p = conn->getPartition();
		delete conn;
	}

	int largest_size = 0;
	std::unordered_set<custom_node> nodes_largest_component;
	for(auto& set_of_partition:p.getSubsets()){
		if(set_of_partition.size() > largest_size){
			largest_size = set_of_partition.size();
			nodes_largest_component.clear();
			for(auto& element:set_of_partition){
				nodes_largest_component.insert(element);
			}
		}
	}

	// delete all edges whose endpoints are not in the node set
	int survived_edges = 0;

	source_graph->forEdges([&](custom_node u, custom_node v, custom_weight w) {
		// if one of the endpoints is not in the node set

			if (nodes_largest_component.find(u) == nodes_largest_component.end() || nodes_largest_component.find(v) == nodes_largest_component.end()){
				if(!source_graph->isDirected())
					assert(nodes_largest_component.find(u) == nodes_largest_component.end() && nodes_largest_component.find(v) == nodes_largest_component.end());
				source_graph->removeEdge(u,v);
			}

			else
				survived_edges++;




	});
	// delete all nodes that are not in the node set
	source_graph->forNodes([&](custom_node u) {
		if (nodes_largest_component.find(u) == nodes_largest_component.end()) {
			source_graph->removeNode(u);
		}

	});

	assert(largest_size==source_graph->numberOfNodes());
	assert(survived_edges==source_graph->numberOfEdges());
	assert(source_graph->hasEdgeIds());
	std::cout<<"After SCC, graph has "<<source_graph->numberOfNodes()<<" vertices and "<<source_graph->numberOfEdges()<<" edges"<<std::endl;
	assert(source_graph->hasEdgeIds());
	source_graph->removeSelfLoops();
	assert(source_graph->numberOfSelfLoops()==0);



}

void Auxiliary::real_world_parser(std::string source, bool weighted, bool directed, custom_node num_mod_tot, custom_node mods_selector){

	std::string output = source+".hist";
	std::string output_mods = source+".mods";
	std::fstream ofs;
	std::fstream ofs_mods;
	std::ifstream ifs(source);
	if (!ifs)
			throw std::runtime_error("Error opening File ");
	ofs.open(output,std::ios::out);
	ofs_mods.open(output_mods,std::ios::out);

	std::string line,temp_field;
	std::string commands = " wc -l " + source +" | awk '{print $1}'";

	int n_lines = atoi(exec(commands.c_str()).c_str());


	ProgressStream sorter(n_lines);
	sorter.label() << "Sorting file "<<source<< " with " << n_lines << " lines w.r.t. time";
	getline(ifs,line);//throw first line away
	++sorter;
	custom_node n1, n2;
	std::string tmp_weight;
	custom_time tmp_time;
	int jumps = 0;
	std::vector<std::tuple<custom_time,custom_node,custom_node,std::string>> rows;
	while (true){

		getline(ifs,line);

		if(ifs.eof())
			break;
		++sorter;
		std::stringstream ss(line);
		ss >> n1 >> n2 >> tmp_weight >> tmp_time;
		if(n1==n2){
			++jumps;
			continue;
		}
		rows.push_back(std::make_tuple(tmp_time,n1,n2,tmp_weight));
	}

	std::sort(rows.begin(),rows.end());
	assert(rows.size() == n_lines- 1 - jumps);

	ProgressStream reader(rows.size());
	reader.label() << "Converting " << rows.size() << " tuples to Real_World+Hist format";

	custom_node v1, v2;
	std::string weight;
	custom_time time;

	std::set<custom_node> vertices;
	std::set<std::pair<custom_node,custom_node> > starting_edges;
	std::vector<NetworKit::GraphEvent> operations;
	std::map<std::pair<custom_node,custom_node>, custom_weight> weightsMapping;
	std::map<std::pair<custom_node,custom_node>, custom_time> timesMapping;

	int insertions = 0;
	int deletions = 0;

	for(int i = 0;i<rows.size();i++){


		++reader;
		v1 = std::get<1>(rows[i]);
		v2 = std::get<2>(rows[i]);
		weight = std::get<3>(rows[i]);
		time = std::get<0>(rows[i]);


		if(i >= rows.size() - num_mod_tot){
			NetworKit::GraphEvent event;
			event.u = v1;
			event.v = v2;
			event.w = boost::lexical_cast<custom_weight>(weight);
			if(weight.compare("-1") == 0 && (mods_selector==0 || mods_selector==2)){
				event.type =  NetworKit::GraphEvent::EDGE_REMOVAL;
				operations.push_back(event);deletions++;
			}
			else if((weight.compare("+1") == 0 || weight.compare("1") == 0) && (mods_selector==1 || mods_selector==2)){
				event.type =  NetworKit::GraphEvent::EDGE_ADDITION;
				operations.push_back(event);insertions++;
			}
			else{
				continue;
			}
			vertices.insert(v1);
			vertices.insert(v2);
		}
		else{
			if(!directed)
				if(starting_edges.find(std::make_pair(v1,v2))!=starting_edges.end() || starting_edges.find(std::make_pair(v2,v1))!=starting_edges.end())
					continue;

			if(directed)
				if(starting_edges.find(std::make_pair(v1,v2))!=starting_edges.end())
					continue;

			vertices.insert(v1);
			vertices.insert(v2);
			starting_edges.insert(std::make_pair(v1,v2));

			weightsMapping[std::make_pair(v1,v2)] = boost::lexical_cast<custom_weight>(weight);
			timesMapping[std::make_pair(v1,v2)] = time;
		}
	}

	std::map<custom_node,custom_node> remapping;

	//remapping vertices to 0...n scale

	custom_node count = 0;
	for(std::set<custom_node>::iterator iter = vertices.begin();iter!=vertices.end();iter++){
	  remapping[*iter] = count;
	  count++;
	}
	if(count!=vertices.size())
	  throw std::runtime_error("wrong counting");


	ofs<<vertices.size()<<" "<<starting_edges.size()<<" "<<weighted<<" "<<directed<<"\n";
	ProgressStream writer(starting_edges.size());
	writer.label() << "Writing " << starting_edges.size() << " edges ";
	for(std::set<std::pair<custom_node,custom_node> >::iterator iter = starting_edges.begin(); iter!=starting_edges.end(); iter++){
		ofs << boost::lexical_cast<std::string>(timesMapping[std::make_pair((*iter).first,(*iter).second)]) << " "
			<< boost::lexical_cast<std::string>(remapping[(*iter).first]) << " "
			<< boost::lexical_cast<std::string>(remapping[(*iter).second]) << " "
			<< boost::lexical_cast<std::string>(weightsMapping[std::make_pair((*iter).first,(*iter).second)])<<std::endl;
		++writer;
	}
	assert(operations.size()==num_mod_tot);
	assert(deletions+insertions==operations.size());
	for(auto& el : operations){
		el.u = remapping[el.u];
		el.v = remapping[el.v];
	}

	ofs_mods<<insertions<<" "<<deletions<<std::endl;
	for(auto& el : operations){
		if(el.type==NetworKit::GraphEvent::EDGE_ADDITION)
			ofs_mods << "+ "<< boost::lexical_cast<std::string>(el.u) <<" "<<  boost::lexical_cast<std::string>(el.v)  <<" "<<boost::lexical_cast<std::string>(el.w)<<std::endl;
		else
			ofs_mods << "- "<< boost::lexical_cast<std::string>(el.u) <<" "<<  boost::lexical_cast<std::string>(el.v)  <<" "<<boost::lexical_cast<std::string>(el.w)<<std::endl;
	}




	ofs.close();
	ofs_mods.close();
	std::cout<<"Graph has initially "<<starting_edges.size()<<" edges, "<<deletions<<" delete and "<<insertions<<" insert operations can be performed"<<std::endl;

}


void Auxiliary::convert(std::string source_file_path, unsigned int type, unsigned int preamble_length, bool weighted, bool directed){

	 std::string output = source_file_path+".hist";
	 std::fstream ofs;
	 ofs.open(output,std::ios::out);
	 std::string line,temp_field;
	 std::string commands = " wc -l " + source_file_path +" | awk '{print $1}'";

	 int n_lines = atoi(exec(commands.c_str()).c_str());
	 std::ifstream ifs(source_file_path);
	 if (!ifs)
	 	throw std::runtime_error("Error opening File ");

	 //"Type of Source Graph [SIGNED(0) UNSIGNED(1) BGR(2) DISTRIBUTED(3) KONECT(4) SNAP(5)]"


	 if(type == 0 && weighted)
		 throw std::runtime_error("Asking to read a signed graph with weights");
	 if(type == 1 && weighted)
		 throw std::runtime_error("Asking to read an unsigned graph with weights");

	 switch(type){
	 	 case 0:{
	 		 custom_node v1, v2, ignored_int;
	 		 ProgressStream reader(n_lines);
	 		 reader.label() << "Converting file with " << n_lines << " lines from Signed to Hist format";

	 		 unsigned int skipped = 0;
	 		 while(skipped<preamble_length){
				  getline(ifs,line);++reader;
				  skipped++;
	 		 }
	 		 std::set<custom_node> vertices;
	 		 std::set<std::pair<custom_node,custom_node> > edges;

	 		 while (true){
	 			 getline(ifs,line);

				  if(ifs.eof())
					break;
				  ++reader;
				  std::stringstream ss(line);
				  ss >> v1 >> v2 >> ignored_int;//can we use the sign?

				  //remove fake edges
				  if(v1==v2)
					continue;
				  if(!directed){
					  if(edges.find(std::make_pair(v1,v2))!=edges.end() || edges.find(std::make_pair(v2,v1))!=edges.end())
						  continue;
				  }
				  if(directed){
					  if(edges.find(std::make_pair(v1,v2))!=edges.end())
						  continue;
				  }
				  vertices.insert(v1);
				  vertices.insert(v2);
				  edges.insert(std::make_pair(v1,v2));
	 		 }

	 		 std::map<custom_node,custom_node> remapping;
	 		 //remapping vertices to 0...n scale
	 		 custom_node count = 0;
	 		 for(std::set<custom_node>::iterator iter = vertices.begin();iter!=vertices.end();iter++){
	 			 remapping[*iter] = count;
	 			 count++;
	 		 }
	 		 if(count!=vertices.size())
	 			 throw std::runtime_error("wrong counting");

	 		 ofs<<vertices.size()<<" "<<edges.size()<<" "<<weighted<<" "<<directed<<"\n";
	 		 for(std::set<std::pair<custom_node,custom_node> >::iterator iter = edges.begin(); iter!=edges.end(); iter++)
	 			 ofs << boost::lexical_cast<std::string>(0) << " "<< boost::lexical_cast<std::string>(remapping[(*iter).first]) <<" "<<  boost::lexical_cast<std::string>(remapping[(*iter).second])  <<" "<<1<<std::endl;

	 		 break;
	 	 }
	 	 case 1:{
	 		 ProgressStream reader(n_lines);
	 		 reader.label() << "Converting file with " << n_lines << " lines from Unsigned (unweighted) to Hist format";


	 		custom_node v1, v2;
			unsigned int skipped = 0;
			while(skipped<preamble_length){
			  getline(ifs,line);++reader;
			  skipped++;
			}
			std::set<custom_node> vertices;
			std::set<std::pair<custom_node,custom_node> > edges;

			while (true){
			  getline(ifs,line);
			  if(ifs.eof())
				break;
			  ++reader;
			  std::stringstream ss(line);
			  ss >> v1 >> v2;


			  //remove fake edges
			  if(v1==v2)
				continue;
			  if(!directed){
				  if(edges.find(std::make_pair(v1,v2))!=edges.end() || edges.find(std::make_pair(v2,v1))!=edges.end())
					  continue;
			  }
			  if(directed){
				  if(edges.find(std::make_pair(v1,v2))!=edges.end())
					  continue;
			  }
			  vertices.insert(v1);
			  vertices.insert(v2);
			  edges.insert(std::make_pair(v1,v2));
			}
			std::map<custom_node,custom_node> remapping;
			//remapping vertices to 0...n scale
			custom_node count = 0;
			for(std::set<custom_node>::iterator iter = vertices.begin();iter!=vertices.end();iter++){
			  remapping[*iter] = count;
			  count++;
			}
			if(count!=vertices.size())
			  throw std::runtime_error("wrong counting");

			ofs<<vertices.size()<<" "<<edges.size()<<" "<<weighted<<" "<<directed<<"\n";
			for(std::set<std::pair<custom_node,custom_node> >::iterator iter = edges.begin(); iter!=edges.end(); iter++)
			  ofs << boost::lexical_cast<std::string>(0) << " "<< boost::lexical_cast<std::string>(remapping[(*iter).first]) <<" "<<  boost::lexical_cast<std::string>(remapping[(*iter).second])  <<" "<<1<<std::endl;

			break;
	 	 }
	 	 case 2:{

	 		 int bgr_type;
	 		 readPrimitive(ifs, bgr_type);
	 		 if (bgr_type != 1)
	 			 throw std::runtime_error("Error while reading BGR File");

	 		 unsigned int readNodes,readEdges;
	 		 readPrimitive(ifs, readNodes);
	 		 readPrimitive(ifs, readEdges);
	 		 std::set<custom_node> vertices;
	 		 std::set<std::pair<custom_node,custom_node> > edges;
	 		 std::map<std::pair<custom_node,custom_node>, custom_weight> weightsMapping;
	 		 ProgressStream reader(readNodes+readEdges);
	 		 reader.label() << "Converting file containing " <<readNodes<<" nodes and "<<readEdges << " edges from BGR to Hist format";


	 		 int coord;
	 		 //ignore coordinates
	 		 for (unsigned int i=0; i < readNodes; ++i){
				readPrimitive(ifs, coord);
				readPrimitive(ifs, coord);
				++reader;
	 		 }

	 		 for (unsigned int i=0; i < readEdges; ++i){
			  DEDGE edge;
			  readPrimitive(ifs, edge);
			  ++reader;



			  //remove fake edges
			  if(edge.source==edge.target)
				continue;

			  if(!directed){
				  if(edges.find(std::make_pair((custom_node)edge.source,(custom_node)edge.target))!=edges.end() || edges.find(std::make_pair((custom_node)edge.target,(custom_node)edge.source))!=edges.end())
					  continue;
			  }
			  if(directed){
				  if(edges.find(std::make_pair((custom_node)edge.source,(custom_node)edge.target))!=edges.end())
					  continue;
			  }
			  vertices.insert((custom_node)edge.source);
			  vertices.insert((custom_node)edge.target);
			  edges.insert(std::make_pair((custom_node)edge.source,(custom_node)edge.target));
			  weightsMapping[std::make_pair((custom_node)edge.source,(custom_node)edge.target)] = (custom_weight)edge.weight;

			}
			std::map<custom_node,custom_node> remapping;
			//remapping vertices to 0...n scale
			custom_node count = 0;
			for(std::set<custom_node>::iterator iter = vertices.begin();iter!=vertices.end();iter++){
			  remapping[*iter] = count;
			  count++;
			}
			if(count!=vertices.size()){
			  throw std::runtime_error("wrong counting");
			}

			ofs<<vertices.size()<<" "<<edges.size()<<" "<<weighted<<" "<<directed<<"\n";
			for(std::set<std::pair<custom_node,custom_node> >::iterator iter = edges.begin(); iter!=edges.end(); iter++)
			  ofs << boost::lexical_cast<std::string>(0) << " "<< boost::lexical_cast<std::string>(remapping[(*iter).first]) <<" "<<  boost::lexical_cast<std::string>(remapping[(*iter).second])  <<" "<<boost::lexical_cast<std::string>(weightsMapping[std::make_pair((*iter).first,(*iter).second)])<<std::endl;

			break;
	 	 }
	 	 case 3:{
			ProgressStream reader(n_lines);
			reader.label() << "Converting file with " << n_lines << " lines from Distributed to Hist format";
			custom_node v1, v2;


			getline(ifs,line);++reader;

			std::set<custom_node> vertices;
			std::set<std::pair<custom_node,custom_node> > edges;
			std::map<std::pair<custom_node,custom_node>, custom_weight> weightsMapping;

			custom_weight we;
			while (true){
			  getline(ifs,line);

			  if(ifs.eof())
				break;
			  ++reader;
			  std::stringstream ss(line);
			  ss >> v1 >> v2 >> we >> we;//can we use the sign?

			  //remove fake edges
			  if(v1==v2)
				continue;

			  if(!directed)
				  if(edges.find(std::make_pair(v1,v2))!=edges.end() || edges.find(std::make_pair(v2,v1))!=edges.end())
					continue;

			  if(directed)
				  if(edges.find(std::make_pair(v1,v2))!=edges.end())
					continue;

			  vertices.insert(v1);
			  vertices.insert(v2);
			  edges.insert(std::make_pair(v1,v2));
			  weightsMapping[std::make_pair(v1,v2)] = we;

			}

			std::map<custom_node,custom_node> remapping;
			//remapping vertices to 0...n scale
			custom_node count = 0;
			for(std::set<custom_node>::iterator iter = vertices.begin();iter!=vertices.end();iter++){
			  remapping[*iter] = count;
			  count++;
			}
			if(count!=vertices.size())
			  throw std::runtime_error("wrong counting");

			ofs<<vertices.size()<<" "<<edges.size()<<" "<<weighted<<" "<<directed<<"\n";
			for(std::set<std::pair<custom_node,custom_node> >::iterator iter = edges.begin(); iter!=edges.end(); iter++)
			  ofs << boost::lexical_cast<std::string>(0) << " "<< boost::lexical_cast<std::string>(remapping[(*iter).first]) <<" "<<  boost::lexical_cast<std::string>(remapping[(*iter).second])  <<" "<<boost::lexical_cast<std::string>(weightsMapping[std::make_pair((*iter).first,(*iter).second)])<<std::endl;
			break;
	 	 }


	 	 case 4:{
	 		std::string formatString = " KONECT";

	 		NetworKit::Graph* graph = new NetworKit::Graph();
			std::cout<<"Reading graph stored in " << source_file_path <<formatString<< " FORMAT"<<std::endl;
			NetworKit::KONECTGraphReader* reader = new NetworKit::KONECTGraphReader();
			*graph = reader->read(source_file_path);
			std::string t1 = graph->isWeighted() ? "weighted " : "unweighted ";
			std::string t2 = graph->isDirected() ? "directed " : "undirected ";
			std::cout<<"Found " <<  t1 << t2  <<"graph in "<<source_file_path<<" containing " << graph->numberOfNodes() << " vertices and " << graph->numberOfEdges() << " edges\n";

			ofs<<graph->numberOfNodes()<<" "<<graph->numberOfEdges()<<" "<<graph->isWeighted()<<" "<<graph->isDirected()<<"\n";
			graph->forEdges([&](custom_node u, custom_node v, custom_weight w) {
				ofs << boost::lexical_cast<std::string>(0) << " "<< boost::lexical_cast<std::string>(u) <<" "<<  boost::lexical_cast<std::string>(v)  <<" "<<boost::lexical_cast<std::string>(w)<<std::endl;
			});

			delete reader;
			delete graph;
			break;
	 	 }

	 	case 5:{
	 		std::string formatString = " SNAP";

			NetworKit::Graph* graph = new NetworKit::Graph();
			std::cout<<"Reading graph stored in " << source_file_path <<formatString<< " FORMAT"<<std::endl;
			NetworKit::SNAPGraphReader* reader = new NetworKit::SNAPGraphReader();
			*graph = reader->read(source_file_path);
			std::string t1 = graph->isWeighted() ? "weighted " : "unweighted ";
			std::string t2 = graph->isDirected() ? "directed " : "undirected ";
			std::cout<<"Found " <<  t1 << t2  <<"graph in "<<source_file_path<<" containing " << graph->numberOfNodes() << " vertices and " << graph->numberOfEdges() << " edges\n";

			ofs<<graph->numberOfNodes()<<" "<<graph->numberOfEdges()<<" "<<graph->isWeighted()<<" "<<graph->isDirected()<<"\n";
			graph->forEdges([&](custom_node u, custom_node v, custom_weight w) {
				ofs << boost::lexical_cast<std::string>(0) << " "<< boost::lexical_cast<std::string>(u) <<" "<<  boost::lexical_cast<std::string>(v)  <<" "<<boost::lexical_cast<std::string>(w)<<std::endl;
			});

			delete reader;
			delete graph;
			break;
		 }


	 	 default:{
	 		 throw std::runtime_error("unknown type_of_source_file");
	 		 break;
	 	 }

	 }
	 ofs.close();
}


double Auxiliary::avg(std::vector<double>& v){
  double value = 0;
  for(size_t t = 0; t < v.size();t++){
    value+=v[t];
  }
  value/=v.size();
  return value;
}
double Auxiliary::avg(std::vector<int>& v){
  double value = 0;
  for(size_t t = 0; t < v.size();t++){
    value+=v[t];
  }
  value/=v.size();
  return value;
}

double Auxiliary::stddev(std::vector<double>& v){
  double average = avg(v);
  double variation = 0;
  for(size_t t = 0; t < v.size();t++)
    variation+=((v[t]-average)*(v[t]-average));

  variation/=v.size();
  return sqrt(variation);

}
double Auxiliary::stddev(std::vector<int>& v){
  double average = avg(v);
  double variation = 0;
  for(size_t t = 0; t < v.size();t++)
    variation+=((v[t]-average)*(v[t]-average));

  variation/=v.size();
  return sqrt(variation);

}

std::string Auxiliary::exec(const char* cmd) {
    std::shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
    if (!pipe)
    	return "ERROR";

    char buffer[128];
    std::string result = "";

    while (!feof(pipe.get())) {
        if (fgets(buffer, 128, pipe.get()) != NULL)
            result += buffer;
    }
    return result;
}

bool Auxiliary::exists (const std::string& name) {
	std::ifstream f(name.c_str());
	return f.good();
}

void Auxiliary::generate(NetworKit::Graph** g,int bound_size, bool directed){
	MersenneTwister* randomGenerator;
	randomGenerator = new MersenneTwister();
	int nodes = round(sqrt(bound_size));
	if(directed){
		INFO("Generating Directed Erdos-Renyi graph");
		NetworKit::ErdosRenyiGenerator* e = new NetworKit::ErdosRenyiGenerator(nodes,(double)(randomGenerator->getRandomInteger()%100+1)/100,1);
		*g = new NetworKit::Graph(e->generate());
		delete e;
		delete randomGenerator;
		return;
	}
	else{
		INFO("Generating Undirected Barabasi-Albert graph");
		NetworKit::BarabasiAlbertGenerator* e = new NetworKit::BarabasiAlbertGenerator(2,nodes,2);
		*g = new NetworKit::Graph(e->generate());
		delete e;
		delete randomGenerator;
		return;

	}




}

void Auxiliary::store_graph(std::string source, NetworKit::Graph** g){

	std::ofstream ofs(source);

	if (!ofs)
		throw std::runtime_error("Error opening File ");

	ofs << (*g)->numberOfNodes() << " " << (*g)->numberOfEdges() << " "
			<<  (*g)->isWeighted() << " " <<  (*g)->isDirected()<<"\n";

	ProgressStream writer((*g)->numberOfEdges());
	std::string t1 =  (*g)->isWeighted()==0 ? " unweighted" : " weighted";
	std::string t2 = (*g)->isDirected()==0 ? " undirected" : " directed";
	writer.label() << "Writing "<< t1 << t2 << " graph in "<<source<<" (CUSTOM FORMAT) containing " <<  (*g)->numberOfNodes()  << " vertices and " << (*g)->numberOfEdges() << " edges ";

	(*g)->forEdges([&](custom_node u, custom_node v, custom_weight w) {
		ofs<<0<<" "<<u<<" "<<v<<" "<<w<<"\n";
		++writer;
	});

	ofs.close();





}
void Auxiliary::copy(NetworKit::Graph** source,NetworKit::Graph** target){
	*target = new NetworKit::Graph((*source)->numberOfNodes(),(*source)->isWeighted(),(*source)->isDirected());
	(*source)->forEdges([&](custom_node u, custom_node v, custom_weight w) {
		(*target)->addEdge(u,v,w);
	});



}

void Auxiliary::read(std::string source, bool scc_, NetworKit::Graph** g){


	std::ifstream ifs(source);

	if (!ifs)
		throw std::runtime_error("Error opening File ");




	int vertices = -1, edges = -1, weighted = -1, directed = -1;
	ifs >> vertices >> edges >> weighted >> directed;
	assert((weighted==0||weighted==1) && (directed==0||directed==1) && (vertices>=0 && edges>=0));
	ProgressStream reader(edges);
	std::string t1 = weighted==0 ? " unweighted" : " weighted";
	std::string t2 = directed==0 ? " undirected" : " directed";
	reader.label() << "Reading"<< t1 << t2 << " graph in "<<source<<" (CUSTOM FORMAT) containing " << vertices << " vertices and " << edges << " edges ";
	NetworKit::Graph* graph = new NetworKit::Graph(vertices,weighted,directed);
	int time,v1,v2,weight;

	while(true){

		ifs >> time >> v1 >> v2 >> weight;
		if(ifs.eof())
			break;

		++reader;

		assert(weighted || weight == 1 || weight == -1);

		if(v1==v2)
			continue;

		assert(graph->hasNode(v1) && graph->hasNode(v2));
		if(graph->hasEdge(v1,v2))
			std::cout<<"SKIPPING ROW"<<std::endl;
		else{
			graph->addEdge(v1,v2,weight);
			#ifndef NDEBUG
				if(!directed){
					if(!graph->hasEdge(v1,v2) && !graph->hasEdge(v2,v1))
						throw std::runtime_error("wrong edge insertion during construction");
				}
				else{
					if(!graph->hasEdge(v1,v2))
						throw std::runtime_error("wrong edge insertion during construction");
				}
				#endif
		}
	}



	ifs.close();
	graph->indexEdges();
	if(scc_)
		SCC(graph);
	else
		INFO("whole graph considered");
	*g = graph;

}

Auxiliary::~Auxiliary() {
}

std::string Auxiliary::getTimeString() {
	std::time_t rawtime;
	std::tm* timeinfo;
	char buffer [80];

	std::time(&rawtime);
	timeinfo = std::localtime(&rawtime);
	std::strftime(buffer,80,"%d-%m-%Y-%H-%M-%S",timeinfo);
	std::string temp = buffer;
	return temp;
}


