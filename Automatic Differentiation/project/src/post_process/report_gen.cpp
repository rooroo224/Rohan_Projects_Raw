// basic file operations
#include <iostream>
#include <fstream>
//using namespace std;

int prepare_plotter(int num_dim, std::string data_file){

	std::ofstream myfile_out;
	std::ifstream gnu_plot_style; 
	std::ifstream gnu_plotter_line; 

	myfile_out.open ("plt_all");
	gnu_plot_style.open ("./src/tex/gnu_plot_style");
	gnu_plotter_line.open ("./src/tex/gnu_plot_plotter_line");

	std::string line;
	for(int i = 0; i<8; i++){
		getline(gnu_plot_style,line);
		myfile_out << line << "\n";
	}
	myfile_out << "\n";
	
	for(int ii = 0; ii<num_dim; ii++){
  		gnu_plotter_line.seekg (0, std::ios::beg);
  		
		getline(gnu_plotter_line,line);
		line.replace(21, 1, std::to_string(ii)); //png name
		myfile_out << line << "\n";
		
		getline(gnu_plotter_line,line);
		line.replace(26, 1, std::to_string(ii+3));  //data col		
		line.replace(6, 8, data_file);	//data file	
		myfile_out << line << "\n";

		myfile_out << "\n";
	}

	myfile_out.close();
	gnu_plot_style.close();
	gnu_plotter_line.close();

	return 0;
}

int prepare_reporter(int num_dim){

	std::ofstream myfile_out;
	std::ifstream tex_body; 
	std::ifstream tex_inc_plot; 

	myfile_out.open ("new.tex");
	tex_body.open ("./src/tex/tex_body");
	tex_inc_plot.open ("./src/tex/tex_inc_plot");

	std::string line;
	for(int i = 0; i<6; i++){
		getline(tex_body,line);
		myfile_out << line << "\n";
	}
	
	for(int ii = 0; ii<num_dim; ii++){
  		tex_inc_plot.seekg (0, std::ios::beg);
		for(int i = 0; i<2; i++){
			getline(tex_inc_plot,line);
			myfile_out << line << "\n";
		}
		
		getline(tex_inc_plot,line);
		line.replace(49, 1, std::to_string(ii));
		myfile_out << line << "\n";
		
		getline(tex_inc_plot,line);
		line.replace(28, 1, std::to_string(ii));
		myfile_out << line << "\n";
		
		for(int i = 0; i<2; i++){
			getline(tex_inc_plot,line);
			myfile_out << line << "\n";
		}
		myfile_out << "\n";
	}
	
	for(int i = 0; i<2; i++){
		getline(tex_body,line);
		myfile_out << line << "\n";
	}

	myfile_out.close();
	tex_body.close();
	tex_inc_plot.close();
	
	return 0;

}

int report_gen(int num_dim, std::string data_file) {

prepare_plotter(num_dim, data_file);

prepare_reporter(num_dim);

std::cout << "USE: \ngnuplot plt_all && pdflatex new.tex >/dev/null  && rm *.log *.out *.aux" << std::endl;

return 0;
}
