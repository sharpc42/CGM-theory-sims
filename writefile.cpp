/* Code to write a given vector array of string data
   to a file under the given filename, which is then saved.
   User must provide desired extension.
*/

#include <iostream>
#include <fstream>

using namespace std;

/*
string filename;
vector<string> data;
*/

void writefile(string filename, vector<string> data) {

	//filename = filename;
	//data = data;

	cout << "Filename is " << filename;

	ofstream output;
	output.open(filename);

	for (int i = 0; i < data.size(); i++) {
		cout << "Writing " << data[i] << " to file...\n";
 		output << data[i] << endl;
	}

	output.close();

}

int main() {}