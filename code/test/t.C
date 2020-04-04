#include <iostream>
using namespace std;
void t(){
	int len = 2;
	char a[];
	a[0] = 'a';
	a[1] = 'b';
	cout<<a[0]<<endl;
	string b(a,len);
	cout<<b<<endl;
	

}
