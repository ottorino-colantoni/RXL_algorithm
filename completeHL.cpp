//============================================================================
// Name        : completeHL.cpp
// Author      : Mattia D'Emidio
// Version     :
// Copyright   : Copyright © 2016 by Mattia D'Emidio
// Description :
//============================================================================

#include "CustomDataTypes.h"
#include "Auxiliary.h"
#include "Labeling.h"
#include "Labeling_Tools.h"

#include "mytimer.h"
#include "callgrind.h"
#include <unordered_set>


int main(int argc, char** argv) {


	//declare supported options
	namespace po = boost::program_options;
	po::options_description desc("Allowed options");

	desc.add_options()
	("graph_type,a", po::value<int>(), "[RANDOMLY-GEN(0) FROM FILE(1)]")
	("graph_location,g", po::value<std::string>(), "Input Graph File Location (in case of -a 1)")
	("num_mods,k", po::value<int>(), "Number of Updates to Be Performed")
	("mods_source,s", po::value<int>(), "[GENERATED(0) FROM FILE(1)]")
	("mods_location,l", po::value<std::string>(), "Mods File Location (in case of -s 1)")
	("mod_type,m", po::value<int>(), "Type of Mods [DEC(0) INC(1) FUL(2)] (in case of -s 0)")
	("minimality,y",po::value<int>(), "Enable minimality preserving incremental update (in case of -m 1 or -m 2)")
	("num_queries,q", po::value<int>(), "Number of Queries to Be Performed")
	("check_corr,c",po::value<int>(), "Check Correctness [FALSE(0) TRUE(1)]")
	("ordering,o",po::value<int>(), "Type of Node Ordering [DEGREE(0) APPROX-BETWEENESS(1) FILE(2)]")
	("extract_scc,e",po::value<int>(), "Take largest component [FALSE(0) TRUE(1)]")
	("reduce_to_subgraph,r",po::value<int>(), "Take subgraph [FALSE(0) TRUE(1)]")
	("graph_size,b",po::value<int>(), "Bound on the edges of the generated graph/extracted subgraph (in case of -a 0 || -r 1)")
	;


	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify( vm);

	int minimality = -1,graph_type = -1,num_mods = -1, num_queries = -1,  mod_type = -1, check_corr = -1, ordering = -1, extract_scc = -1, mods_source = -1, reduce_to_subgraph = -1, graph_size = -1;
	std::string graph_location = "", mods_location = "";

	if (vm.empty()){
		std::cout << desc << "\n";
		throw std::runtime_error("Empty options");
	}
	if (vm.count("graph_type"))
		graph_type = vm["graph_type"].as<int>();

	if (vm.count("graph_location"))
		graph_location = vm["graph_location"].as<std::string>();

	if (vm.count("mods_location"))
		mods_location = vm["mods_location"].as<std::string>();

	if (vm.count("num_mods"))
		num_mods = vm["num_mods"].as<int>();

	if (vm.count("mod_type"))
		mod_type = vm["mod_type"].as<int>();

	if (vm.count("minimality"))
		minimality = vm["minimality"].as<int>();

	if (vm.count("mods_source"))
		mods_source = vm["mods_source"].as<int>();

	if (vm.count("num_queries"))
		num_queries = vm["num_queries"].as<int>();

	if (vm.count("check_corr"))
		check_corr = vm["check_corr"].as<int>();

	if (vm.count("ordering"))
		ordering = vm["ordering"].as<int>();


	if (vm.count("extract_scc"))
		extract_scc = vm["extract_scc"].as<int>();

	if (vm.count("reduce_to_subgraph"))
		reduce_to_subgraph = vm["reduce_to_subgraph"].as<int>();

	if (vm.count("graph_size"))
		graph_size = vm["graph_size"].as<int>();


	if(graph_type!=0 && graph_type!=1){
			std::cout << desc << "\n";
			throw std::runtime_error("wrong graph_type");
	}


	if(graph_location == "" && graph_type==1){
		std::cout << desc << "\n";
	    throw std::runtime_error("wrong graph_location");
	}
	if(minimality != 0 && minimality != 1){
		std::cout << desc << "\n";
		throw std::runtime_error("wrong minimality");
	}

	if(ordering != 0 && ordering != 1 && ordering != 2){
		std::cout << desc << "\n";
		throw std::runtime_error("wrong ordering");
	}



	if(num_mods<1){
	  std::cout << desc << "\n";
	  throw std::runtime_error("wrong num_mods");
	}

	if(num_queries<1){
	  std::cout << desc << "\n";
	  throw std::runtime_error("wrong num_queries");
	}

	if(mods_source!=0 && mods_source!=1){
	  std::cout << desc << "\n";
	  throw std::runtime_error("wrong mods_source");
	}

	if(mod_type!=0 && mod_type!=1 && mod_type!=2 && mods_source==0){
	  std::cout << desc << "\n";
	  throw std::runtime_error("wrong mod_type");
	}


	if(reduce_to_subgraph!=0 && reduce_to_subgraph!=1){
		std::cout << desc << "\n";
		throw std::runtime_error("wrong reduce_to_subgraph");
	}
	if((reduce_to_subgraph==1 || graph_type==0) && graph_size<=2 ){
		std::cout << desc << "\n";
		throw std::runtime_error("wrong graph_size");
	}


	if(mods_location == "" && mods_source==1){
		std::cout << desc << "\n";
		throw std::runtime_error("wrong mods_location");
	}


	if(check_corr!=0 && check_corr!=1){
		std::cout << desc << "\n";
		throw std::runtime_error("wrong check_corr");
	}

	if(extract_scc!=0 && extract_scc!=1){
	  std::cout << desc << "\n";
	  throw std::runtime_error("wrong extract_scc");
	}



	
	std::vector<double> num_affected,num_affected_x,num_affected_y;
	std::vector<double> num_added_labels,num_removed_labels;
	std::vector<double> V,E;
	std::vector<double> fs_prep_times,dyn_upd_times,ratio_times;
	std::vector<double> fs_sizes,dyn_sizes,ratio_sizes;
	std::vector<double> fs_query_times, dyn_query_times, ratio_query_times;

	std::string workload_string;

	switch(mod_type){
	case (0):{
		if(mods_source==1)
			workload_string="D_FIL";
		else
			workload_string = "DEC";
		break;
	}

	case (1):{
		if(mods_source==1)
			workload_string="I_FIL";
		else
			workload_string = "INC";
		break;
	}
	case (2):{
		if(mods_source==1)
			workload_string="F_FIL";
		else
			workload_string = "FUL";
		break;
	}

	default:
		break;
	}

	std::string order_string;

	switch(ordering){
			case (0):
				order_string = "DEG";
				break;

			case (1):
				order_string = "BET";
				break;
			case (2):
				order_string = "MAN";
				break;
			default:
			break;
	}

	std::string shortenedName = graph_type==1 ? graph_location.substr(0, 16) : "randomly";

	std::string timestampstring = shortenedName+"_"+workload_string+"_"+order_string+"_"+Auxiliary::getTimeString();
	std::string logFile = timestampstring +".res";
	std::ofstream ofs;
	ofs.open(logFile.c_str(),std::ofstream::out);
	ofs<<"#;|V|;|E|;FST;DynT;SpU;FSQT(ms);DynQT(ms);RatQT;FSSize;DynSize;RSize;|aff|;|aff_x|;|aff_y|;Added;Removed"<<std::endl;

	#ifdef _OPENMP
		std::cout<<"Available Threads: "<<omp_get_max_threads()<<std::endl;
	#else
		std::cout<<"Available Threads: "<<1<<std::endl;
	#endif


	NetworKit::Graph* graph;
	if(graph_type==1)
		Auxiliary::read(graph_location,extract_scc,&graph);
	else{
		MersenneTwister* randomGenerator;
		randomGenerator = new MersenneTwister();
		Auxiliary::generate(&graph,graph_size,randomGenerator->getRandomInteger()%2==0);
		std::string rand_name = shortenedName + "_"+Auxiliary::getTimeString()+".hist";
		int c = 1;
		while(Auxiliary::exists(rand_name)){
			rand_name = shortenedName + "_"+Auxiliary::getTimeString()+"_"+boost::lexical_cast<std::string>(c)+".hist";
			++c;
		}

		Auxiliary::store_graph(rand_name,&graph);
		delete randomGenerator;


	}


	NetworKit::Graph* mirror_graph;

	if(check_corr){
		if(graph_type==1)
			Auxiliary::read(graph_location,extract_scc,&mirror_graph);
		else{
			Auxiliary::copy(&graph,&mirror_graph);

		}
		assert(graph->size()==mirror_graph->size());
	}

	if(graph_type==0 && reduce_to_subgraph==1)
		throw new std::runtime_error("asking subgraph of generated graph");

	if(reduce_to_subgraph==1){
		if(graph_size>graph->numberOfEdges())
			throw new std::runtime_error("too large subgraph");
		INFO("Extracting desired subgraph");
		int n_ = round(100*graph_size);
		std::cout<<"Temptative n_nodes: "<<n_<<"\n";
		while(true){
			std::unordered_set<custom_node> nodes;
			for(size_t t=0;t<n_;t++)
				nodes.insert(graph->randomNode());
			NetworKit::Graph red_graph = graph->subgraphFromNodes(nodes);
			if(red_graph.numberOfEdges()>3 && red_graph.numberOfEdges()<graph_size){
				*graph = red_graph;
				break;
			}
			else
				n_--;

		}
		std::cout<<"Done!\n";
		std::cout<<graph->numberOfNodes()<<" "<<graph->numberOfEdges()<<" "<<graph->isWeighted()<<" "<<graph->isDirected()<<"\n";

		graph->forEdges([&] (custom_node u, custom_node v){
			std::cout<<u<<" "<<v<<"\n";
		});
		if(check_corr){
			*mirror_graph=*graph;
			assert(graph->size()==mirror_graph->size());
		}
	}

	std::pair<std::vector<custom_node>,std::vector<custom_node>> order_keeper = Auxiliary::compute_ordering(graph, ordering);


	if(mod_type == 0 && num_mods>=graph->numberOfEdges())
		throw new std::runtime_error("Graph has not enough edges");


	std::vector<NetworKit::GraphEvent> mods;

	if(mods_source==0){
		Auxiliary::init_mods(mods, mod_type, num_mods, graph,0.25);
		std::string modsFilen = timestampstring +".mods";
		Auxiliary::store_mods(modsFilen,mods);

	}
	else
		Auxiliary::read_mods(mods_location, mods);

	assert(mods.size()==num_mods);

	Labeling* labeling = new Labeling(graph->isDirected());
	Labeling_Tools* manager = new Labeling_Tools(graph,labeling,order_keeper);

#ifndef NDEBUG
//	Verify build phase
	Auxiliary::test_match(graph,labeling,100);
#endif

	double fs_size_pre_corr = labeling->getNumberOfLabelEntries();
	double fs_time_pre_corr = manager->preprocessing_time;
	double fs_qtime_pre_corr = 0.0;

	if(!check_corr){
		mytimer t_counter;

		std::cout<< "Measuring query time of from scratch with " <<num_queries<< " correctness queries";
#pragma omp parallel for schedule(guided, 1)
		for(int i = 0;i<num_queries;i++){
			custom_node u = graph->randomNode();
			custom_node v = graph->randomNode();

			t_counter.restart();
			custom_weight dfs = labeling->query(u,v);
			fs_qtime_pre_corr += t_counter.elapsed()*10e3;

		}
#pragma omp flush

		fs_qtime_pre_corr/=num_queries;
		std::cout<<"... done!\n";
	}


	for (size_t t = 0; t<mods.size();++t) {

		NetworKit::GraphEvent event = mods[t];
		assert(event.w!=0);

		if(event.type==NetworKit::GraphEvent::EDGE_REMOVAL){ // removal
			assert(event.w==-1);
			std::cout<<"Performing update " << t+1 << " of "<< num_mods <<" -- Removing edge("<<event.u<<","<<event.v<<")"<<std::endl;
		}
		else if(event.type==NetworKit::GraphEvent::EDGE_ADDITION){

			std::cout<<"Performing update " << t+1 << " of "<< num_mods <<" -- Adding edge("<<event.u<<","<<event.v<<") New_Weight: "<<event.w<<std::endl;
		}
		else if(event.type==NetworKit::GraphEvent::EDGE_WEIGHT_UPDATE && event.w>0){
#ifndef NDEBUG
			if(event.w<=0)
				assert(false);
			if(event.w>2*graph->weight(event.u,event.v))
				assert(false);
#endif
			std::cout<<"Performing update " << t+1 << " of "<< num_mods <<" -- Decremental update on edge("<<event.u<<","<<event.v<<")"
									<<" Old_Weight: "<<graph->weight(event.u,event.v)<<" New_Weight: "<<event.w<<std::endl;

		}
		else{
			assert(event.w<0);
			assert(event.type==NetworKit::GraphEvent::EDGE_WEIGHT_UPDATE);
			assert(-event.w>0);
			assert(graph->weight(event.u,event.v)>-event.w);
			std::cout<<"Performing update " << t+1 << " of "<< num_mods <<" -- Incremental update on edge("<<event.u<<","<<event.v<<")"
					<<" Old_Weight: "<<graph->weight(event.u,event.v)<<" New_Weight: "<<-event.w<<std::endl;
		}

		UpdateData* upd_d = new UpdateData();

		manager->update(&event,upd_d,minimality); //opera su una struct con i dati dell'aggiornamento

		double dyn_time = upd_d->time;
		double dyn_size = labeling->getNumberOfLabelEntries();

		dyn_upd_times.push_back(dyn_time);
		dyn_sizes.push_back(dyn_size);
		num_affected.push_back(upd_d->affected);
		num_affected_x.push_back(upd_d->affected_x);
		num_affected_y.push_back(upd_d->affected_y);
		num_added_labels.push_back(upd_d->added);
		num_removed_labels.push_back(upd_d->removed);
		V.push_back(graph->numberOfNodes());
		E.push_back(graph->numberOfEdges());

		delete upd_d;


		if(check_corr){
				std::vector<NetworKit::GraphEvent> temp;

				NetworKit::GraphUpdater graph_updater(*mirror_graph);
				temp.push_back(event);
				if(NetworKit::GraphEvent::EDGE_WEIGHT_UPDATE && event.w<0)
					temp[0].w = -temp[0].w;//incremental wchange case

				assert(temp[0].w>0);

				graph_updater.update(temp);

				temp.clear();

				Labeling* labeling_fromscratch = new Labeling(mirror_graph->isDirected());
				Labeling_Tools* manager_fromscratch = new Labeling_Tools(mirror_graph,labeling_fromscratch,order_keeper,false);

				assert(graph->numberOfNodes()==mirror_graph->numberOfNodes());
				assert(graph->numberOfEdges()==mirror_graph->numberOfEdges());


				fs_prep_times.push_back(manager_fromscratch->preprocessing_time);
				ratio_times.push_back(manager_fromscratch->preprocessing_time/dyn_time);
				fs_sizes.push_back(labeling_fromscratch->getNumberOfLabelEntries());
				ratio_sizes.push_back(labeling_fromscratch->getNumberOfLabelEntries()/dyn_size);

				double avg_dyn_qtime = 0.0;
				double avg_fs_qtime = 0.0;


				std::cout<< "Performing " <<num_queries<< " correctness queries";
				mytimer time_counter;



				for(int i = 0;i<num_queries;i++){
					custom_node u = graph->randomNode();
					custom_node v = graph->randomNode();

					custom_weight d_dyn,d_fs;
					if(i%2==0){
						time_counter.restart();
						d_dyn = labeling->query(u,v);
						avg_dyn_qtime += time_counter.elapsed()*10e3;
						time_counter.restart();
						d_fs = labeling_fromscratch->query(u,v);
						avg_fs_qtime += time_counter.elapsed()*10e3;
					}
					else{
						time_counter.restart();
						d_fs = labeling_fromscratch->query(u,v);
						avg_fs_qtime += time_counter.elapsed()*10e3;
						time_counter.restart();
						d_dyn = labeling->query(u,v);
						avg_dyn_qtime += time_counter.elapsed()*10e3;
					}


					if(d_dyn!=d_fs){
						const NetworKit::Graph & G = *graph;
						custom_weight d_graph;
						custom_weight d_mirror_graph;
						if(graph->isWeighted()){
							NetworKit::Dijkstra* dij = new NetworKit::Dijkstra(G,u,true,false,v);
							dij->run();
							d_graph = dij->distance(v);
							const NetworKit::Graph & G_prime = *mirror_graph;
							NetworKit::Dijkstra* dij_prime = new NetworKit::Dijkstra(G_prime,u,true,false,v);
							dij_prime->run();
							d_mirror_graph = dij_prime->distance(v);
						}
						else{
							NetworKit::BFS* dij = new NetworKit::BFS(G,u,true,false,v);
							dij->run();
							d_graph = dij->distance(v);
							const NetworKit::Graph & G_prime = *mirror_graph;
							NetworKit::BFS* dij_prime = new NetworKit::BFS(G_prime,u,true,false,v);
							dij_prime->run();
							d_mirror_graph = dij_prime->distance(v);
						}




						std::cout<<"\nu: "<<u<<" v: "<<v
								<<"\nd_fs: "<<d_fs
								<<"\nd_dyn: "<<d_dyn
								<<"\nd_graph: "<<d_graph
								<<" d_mirror_graph: "<<d_mirror_graph<<std::endl;

						throw std::runtime_error("wrong update by the dynamic algorithm");
					}

				}

				avg_dyn_qtime /= num_queries;
				avg_fs_qtime /= num_queries;

				dyn_query_times.push_back(avg_dyn_qtime);
				fs_query_times.push_back(avg_fs_qtime);
				ratio_query_times.push_back(avg_fs_qtime/avg_dyn_qtime);
				delete labeling_fromscratch;
				delete manager_fromscratch;
				std::cout<<"... done!\n";


		}
		else{

			fs_prep_times.push_back(fs_time_pre_corr);
			ratio_times.push_back(fs_time_pre_corr/dyn_time);
			fs_sizes.push_back(fs_size_pre_corr);
			ratio_sizes.push_back(fs_size_pre_corr/dyn_size);

			std::cout<< "Measuring query times, performing " <<num_queries<< " queries";
			mytimer time_counter;
			double avg_dyn_qtime = 0.0;
#pragma omp parallel for schedule(guided, 1)

			for(int i = 0;i<num_queries;i++){
				custom_node u = graph->randomNode();
				custom_node v = graph->randomNode();
				time_counter.restart();
				custom_weight d_dyn = labeling->query(u,v);
				avg_dyn_qtime += time_counter.elapsed()*10e3;
			}
#pragma omp flush
			avg_dyn_qtime /= num_queries;
			dyn_query_times.push_back(avg_dyn_qtime);
			fs_query_times.push_back(fs_qtime_pre_corr);
			ratio_query_times.push_back(fs_qtime_pre_corr/avg_dyn_qtime);
			std::cout<<"... done!\n";


		}
		ofs<<t<<std::scientific<<std::setprecision(2)
		<<";"<<V[t]<<";"<<E[t]
				<<";"<<fs_prep_times[t]<<";"<<dyn_upd_times[t]<<";"<<std::scientific<<std::setprecision(2)<<ratio_times[t]<<";"
				<<fs_query_times[t]<<";"<<dyn_query_times[t]<<";"<<ratio_query_times[t]<<";"
				<<std::scientific<<std::setprecision(2)<<fs_sizes[t]<<";"<<dyn_sizes[t]<<";"<<ratio_sizes[t]<<";"
				<<num_affected[t]<<";"<<num_affected_x[t]<<";"<<num_affected_y[t]<<";"<<num_added_labels[t]<<";"<<num_removed_labels[t]<<std::endl;

	}


	std::cout<<"AVG;|V|;|E|;FST;DynT;SpU;FSQT(ms);DynQT(ms);RatQT;FSSize;DynSize;RSize;|aff|;|aff_x|;|aff_y|;Added;Removed"<<std::endl;
	std::cout<<"AVG;"<<std::scientific<<std::setprecision(2)<<Auxiliary::avg(V)<<";"<<Auxiliary::avg(E)<<";"<<Auxiliary::avg(fs_prep_times)<<";"<<Auxiliary::avg(dyn_upd_times)<<std::scientific<<std::setprecision(2)<<";"<<Auxiliary::avg(ratio_times)<<";"
		<<Auxiliary::avg(fs_query_times)<<";"<<Auxiliary::avg(dyn_query_times)<<";"<<Auxiliary::avg(ratio_query_times)<<";"
		<<std::scientific<<std::setprecision(2)<<Auxiliary::avg(fs_sizes)<<";"<<Auxiliary::avg(dyn_sizes)<<";"<<Auxiliary::avg(ratio_sizes)<<";"
		<<Auxiliary::avg(num_affected)<<";"<<Auxiliary::avg(num_affected_x)<<";"<<Auxiliary::avg(num_affected_y)<<";"<<Auxiliary::avg(num_added_labels)<<";"<<Auxiliary::avg(num_removed_labels)
	<<std::endl;

	ofs<<"AVG;|V|;|E|;FST;DynT;SpU;FSQT(ms);DynQT(ms);RatQT;FSSize;DynSize;RSize;|aff|;|aff_x|;|aff_y|;Added;Removed"<<std::endl;
	ofs<<"AVG;"<<std::scientific<<std::setprecision(2)<<Auxiliary::avg(V)<<";"<<Auxiliary::avg(E)<<";"<<Auxiliary::avg(fs_prep_times)<<";"<<Auxiliary::avg(dyn_upd_times)<<";"<<std::scientific<<std::setprecision(2)<<Auxiliary::avg(ratio_times)<<";"
		<<Auxiliary::avg(fs_query_times)<<";"<<Auxiliary::avg(dyn_query_times)<<";"<<Auxiliary::avg(ratio_query_times)<<";"
		<<std::scientific<<std::setprecision(2)<<Auxiliary::avg(fs_sizes)<<";"<<Auxiliary::avg(dyn_sizes)<<";"<<Auxiliary::avg(ratio_sizes)<<";"
		<<Auxiliary::avg(num_affected)<<";"<<Auxiliary::avg(num_affected_x)<<";"<<Auxiliary::avg(num_affected_y)<<";"<<Auxiliary::avg(num_added_labels)<<";"<<Auxiliary::avg(num_removed_labels)<<std::endl;

	std::sort(V.begin(),V.end());
	std::sort(E.begin(),E.end());
	std::sort(fs_prep_times.begin(),fs_prep_times.end());
	std::sort(dyn_upd_times.begin(),dyn_upd_times.end());
	std::sort(ratio_times.begin(),ratio_times.end());
	std::sort(fs_query_times.begin(),fs_query_times.end());
	std::sort(dyn_query_times.begin(),dyn_query_times.end());
	std::sort(ratio_query_times.begin(),ratio_query_times.end());
	std::sort(fs_sizes.begin(),fs_sizes.end());
	std::sort(dyn_sizes.begin(),dyn_sizes.end());
	std::sort(ratio_sizes.begin(),ratio_sizes.end());
	std::sort(num_affected.begin(),num_affected.end());
	std::sort(num_affected_x.begin(),num_affected_x.end());
	std::sort(num_affected_y.begin(),num_affected_y.end());
	std::sort(num_added_labels.begin(),num_added_labels.end());
	std::sort(num_removed_labels.begin(),num_removed_labels.end());

	std::cout<<"MED;|V|;|E|;FST;DynT;SpU;FSQT(ms);DynQT(ms);RatQT;FSSize;DynSize;RSize;|aff|;|aff_x|;|aff_y|;Added;Removed"<<std::endl;
	std::cout<<"MED;"<<std::scientific<<std::setprecision(2)<<V[round(V.size()/2)]<<";"<<E[round(E.size()/2)]<<";"<<fs_prep_times[round(fs_prep_times.size()/2)]<<";"<<dyn_upd_times[round(dyn_upd_times.size()/2)]<<std::scientific<<std::setprecision(2)<<";"<<ratio_times[round(ratio_times.size()/2)]<<";"
		<<fs_query_times[round(fs_query_times.size()/2)]<<";"<<dyn_query_times[round(dyn_query_times.size()/2)]<<";"<<ratio_query_times[round(ratio_query_times.size()/2)]<<";"
		<<std::scientific<<std::setprecision(2)<<fs_sizes[round(fs_sizes.size()/2)]<<";"<<dyn_sizes[round(dyn_sizes.size()/2)]<<";"<<ratio_sizes[round(ratio_sizes.size()/2)]<<";"
		<<num_affected[round(num_affected.size()/2)]<<";"<<num_affected_x[round(num_affected_x.size()/2)]<<";"<<num_affected_y[round(num_affected_y.size()/2)]<<";"<<num_added_labels[round(num_added_labels.size()/2)]<<";"<<num_removed_labels[round(num_removed_labels.size()/2)]
	<<std::endl;

	ofs<<"MED;|V|;|E|;FST;DynT;SpU;FSQT(ms);DynQT(ms);RatQT;FSSize;DynSize;RSize;|aff|;|aff_x|;|aff_y|;Added;Removed"<<std::endl;
	ofs<<"MED;"<<std::scientific<<std::setprecision(2)<<V[round(V.size()/2)]<<";"<<E[round(E.size()/2)]<<";"<<fs_prep_times[round(fs_prep_times.size()/2)]<<";"<<dyn_upd_times[round(dyn_upd_times.size()/2)]<<std::scientific<<std::setprecision(2)<<";"<<ratio_times[round(ratio_times.size()/2)]<<";"
			<<fs_query_times[round(fs_query_times.size()/2)]<<";"<<dyn_query_times[round(dyn_query_times.size()/2)]<<";"<<ratio_query_times[round(ratio_query_times.size()/2)]<<";"
			<<std::scientific<<std::setprecision(2)<<fs_sizes[round(fs_sizes.size()/2)]<<";"<<dyn_sizes[round(dyn_sizes.size()/2)]<<";"<<ratio_sizes[round(ratio_sizes.size()/2)]<<";"
			<<num_affected[round(num_affected.size()/2)]<<";"<<num_affected_x[round(num_affected_x.size()/2)]<<";"<<num_affected_y[round(num_affected_y.size()/2)]<<";"<<num_added_labels[round(num_added_labels.size()/2)]<<";"<<num_removed_labels[round(num_removed_labels.size()/2)]
		<<std::endl;

	ofs.close();


	delete graph;
	delete labeling;
	delete manager;

	if(check_corr)
		delete mirror_graph;



	return EXIT_SUCCESS;

	
  
}
