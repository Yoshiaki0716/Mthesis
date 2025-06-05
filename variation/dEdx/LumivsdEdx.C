#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include "utils.h"
#include "/home/yunagai/RootUtils/AtlasStyle.C"
#include "/home/yunagai/RootUtils/AtlasLabels.C"

double getValueFromFile(const std::string &file_path, int row_number, int column_number)
{

	std::ifstream file(file_path);
	row_number--;
	column_number--;

	std::string line;
	double value = 0.0;

	for (int i = 0; i < row_number; ++i)
	{
		if (!std::getline(file, line))
		{
			std::cerr << "Error: Row number exceeds the number of lines in the file." << std::endl;
			return -1;
		}
	}
	if (std::getline(file, line))
	{
		std::istringstream iss(line);
		std::string cell;
		for (int i = 0; i <= column_number; ++i)
		{
			if (!(iss >> cell))
			{
				std::cerr << "Error: Column number exceeds the number of columns in the file." << std::endl;
				return -1;
			}
			if (i == column_number)
			{
				double value;
				std::istringstream(cell) >> value;

				return value;
			}
		}
	}

	file.close();
	return 0;
}

int getTotalCountFromFile(const std::string &file_path, int imode)
{
	std::ifstream file(file_path);

	std::string line;
	int row_count = 0;
	int col_count = 0;

	while (std::getline(file, line))
	{
		row_count++;
		std::istringstream iss(line);
		std::string cell;
		int local_col_count = 0;
		while (iss >> cell)
		{
			local_col_count++;
		}
		col_count = std::max(col_count, local_col_count);
	}

	file.close();

	if (imode == 0)
	{
		return row_count;
	}
	else if (imode == 1)
	{
		return col_count;
	}
	return -1;
}

int LumivsdEdx()
{
	//	SetAtlasStyle();

	std::string directory = "/home/yunagai/Study/data/makePlot/Result/run2/selection_pVtx/";
	std::vector<std::vector<std::string>> filelists(4);
	std::vector<std::string> patterns = {"fitresults_IBL_*.txt", "fitresults_BLY_*.txt", "fitresults_LY1_*.txt", "fitresults_LY2_*.txt"};

	for (size_t i = 0; i < patterns.size(); ++i)
	{
		std::string command = "ls " + directory + patterns[i] + " > filelist_temp.txt";
		system(command.c_str());

		std::ifstream filelist("filelist_temp.txt");
		std::string line;
		while (std::getline(filelist, line))
		{
			std::string filename = line.substr(line.find_last_of("/") + 1);
			filelists[i].push_back(filename);
		}
		filelist.close();
	}

	// Read the luminosity values from lumitable.txt
	std::ifstream lumitable("lumitable.txt");
	std::map<int, double> luminosityMap; // Map to store runnumber-luminosity pairs
	int runnumber;
	double luminosity;
	while (lumitable >> runnumber >> luminosity)
	{
		luminosityMap[runnumber] = luminosity;
	}
	lumitable.close();

	const int lyr_num = filelists.size();
	const int npoint = filelists[0].size();
	const int col_num1 = 5; // MPV
	const int col_num2 = 6; // MPV_error
	int tot_row = getTotalCountFromFile(directory + filelists[0][0], 0);
	std::vector<double> x_value, y_value, y_value_err;
	std::vector<double> filtered_x_value, filtered_y_value, filtered_y_value_err;
	std::vector<std::vector<TGraphErrors *>> tgraph_vector(lyr_num);

	for (size_t h = 0; h < lyr_num; ++h)
	{
		for (int i = 0; i < tot_row; ++i)
		{
			for (int j = 0; j < npoint; ++j)
			{
				// Get x_values from each filenname.
				std::string filename = filelists[h][j];
				std::string runnumber;
				std::regex pattern("fitresults_\\w+_(\\d+)\\.txt");
				std::smatch match;
				std::regex_search(filename, match, pattern);
				runnumber = match[1];
				x_value.push_back(luminosityMap[std::stoi(runnumber)]);

				// Get y_values.
				std::string file_path = directory + filename;
				y_value.push_back(getValueFromFile(file_path, i + 1, col_num1));
				y_value_err.push_back(getValueFromFile(file_path, i + 1, col_num2));

				std::cout << "File name :" << filename << ", x_value = " << x_value[j] << ",  y_value = " << y_value[j] << ",  y_value err= " << y_value_err[j] << std::endl;

				// Use the data as filtered_y_value when the y_value_err is smaller than 1.0
				if (y_value_err[j] <= 0.01)
				{
					filtered_x_value.push_back(luminosityMap[std::stoi(runnumber)]);
					filtered_y_value.push_back(getValueFromFile(file_path, i + 1, col_num1));
					filtered_y_value_err.push_back(getValueFromFile(file_path, i + 1, col_num2));
				}
				else
				{
					std::cout << "Exclued data : i = " << i << ",File name :" << file_path << ", y_err = " << y_value_err[j] << std::endl;
				}
			}

			TGraphErrors *graph = new TGraphErrors(npoint, x_value.data(), y_value.data(), nullptr, y_value_err.data());
			// TGraphErrors *graph = new TGraphErrors(npoint, filtered_x_value.data(), filtered_y_value.data(), nullptr, filtered_y_value_err.data()); //Plot only reliable datas.
			tgraph_vector[h].push_back(graph);

			x_value.clear();
			y_value.clear();
			y_value_err.clear();
			filtered_x_value.clear();
			filtered_y_value.clear();
			filtered_y_value_err.clear();
		}
	}

	TCanvas *c1 = new TCanvas("c1", "Plot", 800, 600);
	std::vector<TLegend *> legends(lyr_num);
	c1->Divide(2, 2);

	for (size_t h = 0; h < lyr_num; ++h)
	{
		c1->cd(h + 1);
		legends[h] = new TLegend(0.7, 0.7, 0.9, 0.9);

		for (int i = 0; i < tot_row; ++i)
		{
			tgraph_vector[0][0]->SetTitle("run2 IBL Integrated luminosity vs dE/dx; Integrated luminosity [fb^{-1}] ; MPV_dE/dx [MeV g^{-1} cm^{2}]");
			tgraph_vector[1][0]->SetTitle("run2 B-layer Integrated luminosity vs dE/dx; Integrated luminosity [fb^{-1}] ; MPV_dE/dx [MeV g^{-1} cm^{2}]");
			tgraph_vector[2][0]->SetTitle("run2 Layer-1 Integrated luminosity vs dE/dx; Integrated luminosity [fb^{-1}] ; MPV_dE/dx [MeV g^{-1} cm^{2}]");
			tgraph_vector[3][0]->SetTitle("run2 Layer-2 Integrated luminosity vs dE/dx; Integrated luminosity [fb^{-1}] ; MPV_dE/dx [MeV g^{-1} cm^{2}]");

			if (i == 0)
			{
				tgraph_vector[h][i]->GetYaxis()->SetTitleOffset(1.5);
				tgraph_vector[h][i]->SetMarkerStyle(21);
				tgraph_vector[h][i]->SetMarkerSize(0.75);
				tgraph_vector[h][i]->SetMarkerColor(4);
				//	tgraph_vector[h][i]->GetXaxis()->SetRangeUser(0, 160);
				//	tgraph_vector[h][i]->GetYaxis()->SetRangeUser(0.5, 1.2);
				tgraph_vector[h][i]->Draw("ape");
			}
			else
			{
				tgraph_vector[h][i]->SetMarkerStyle(21);
				tgraph_vector[h][i]->SetMarkerSize(0.75);
				tgraph_vector[h][i]->SetMarkerColor(5 + 2 * i);
				tgraph_vector[h][i]->Draw("p same");
				tgraph_vector[h][i]->Draw("e same");
			}
			if (h == 1)
			{
				legends[h]->AddEntry(tgraph_vector[1][i], ("Percentage =  " + std::to_string((i + 1) * 0.1 + 0.4)).c_str(), "p");
			}
			else
			{
				legends[h]->AddEntry(tgraph_vector[h][i], ("Percentage =  " + std::to_string((i + 1) * 0.1)).c_str(), "p");
			}
		}

		// Draw lines which mean the beginning of each year.
		std::vector<double> x_lines = {4.2, 4.2 + 38.5, 4.2 + 38.5 + 50.2, 4.2 + 38.5 + 50.2 + 63.3};
		std::vector<std::string> year_labels = {"2016", "2017", "2018"};
		std::vector<TLine *> lines;

		double ymin = tgraph_vector[h][0]->GetYaxis()->GetXmin();
		double ymax = tgraph_vector[h][0]->GetYaxis()->GetXmax();
		double label_y = ymax - 0.2 * (ymax - ymin); //

		for (size_t i = 0; i < x_lines.size(); ++i)
		{
			// Draw vertical lines
			TLine *line = new TLine(x_lines[i], ymin, x_lines[i], ymax);
			line->SetLineStyle(2);
			line->SetLineColor(kGray + 2);
			line->Draw("same");
			lines.push_back(line);

			// Add year labels
			if (i < year_labels.size())
			{
				double label_x;
				if (i == year_labels.size() - 1)
				{
					// For "2018", place the label just to the right of the left line
					label_x = x_lines[i] + 13.0;
				}
				else
				{
					// For other years, place the label in the center
					label_x = (x_lines[i] + x_lines[i + 1]) / 2.0;
				}
				TLatex *latex = new TLatex(label_x, label_y, year_labels[i].c_str());
				latex->SetTextSize(0.04);
				latex->SetTextColor(kGray + 2);
				latex->SetTextAlign(22);
				latex->Draw("same");
			}
		}

		legends[h]->Draw();
		c1->SetLeftMargin(0.15);
		c1->SetBottomMargin(0.15);
		c1->Update();
		c1->Draw();
		ATLASLabel(0.15, 0.15, "Work in progress");
		c1->Print("plot.pdf");
	}

	c1->Update();
	return 0;
}
