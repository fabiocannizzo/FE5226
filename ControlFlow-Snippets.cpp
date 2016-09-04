#include <iostream>

using namespace std;
//
//int main()
//{
//	double sum = 1.0;
//	double lastTerm = 1.0;
//	for (int i = 1; i <= 100; ++i)
//	{
//		lastTerm /= 3.0;
//		sum += lastTerm;
//	}
//
//	cout << sum << endl;
//	return 0;
//}
//
//int main()
//{
//	int x = 1000;
//	while (x >= 0)
//		x -= 13;
//
//	cout << x << endl;
//	return 0;
//}
//
//int main()
//{
//	int x;
//	for (x = 1000; x >= 0; x -= 13)
//		;
//
//	cout << x << endl;
//	return 0;
//}
//
//int main()
//{
//	int cnt = 0;
//
//	bool keepreading = true;
//	while (keepreading)
//	{
//		int x;
//		cout << "enter a number\n";
//		cin >> x;
//		keepreading = (x % 5 != 0);
//		++cnt;
//	}
//
//	cout << cnt << endl;
//	return 0;
//}

//int main()
//{
//	int x;
//	cout << "enter a number\n";
//	cin >> x;
//
//	int cnt = 1;
//
//	while (x % 5 != 0)
//	{
//		++cnt;
//		cout << "enter a number\n";
//		cin >> x;
//	}
//
//	cout << cnt << endl;
//	return 0;
//}

//int main()
//{
//	int x = 1;
//	int cnt = 0;
//
//	while (x % 5 != 0)
//	{
//		++cnt;
//		cout << "enter a number\n";
//		cin >> x;
//	}
//
//	cout << cnt << endl;
//	return 0;
//}

//int main()
//{
//	int cnt = 0;
//	int x;
//	do
//	{
//		++cnt;
//		cout << "enter a number\n";
//		cin >> x;
//	} while (x % 5 != 0);
//
//	cout << cnt << endl;
//	return 0;
//}

//int main()
//{
//	int inRange = 0;
//	for (int i =0; i < 1000000; ++i) {
//		int x = rand();
//		//inRange += ((x >= 500) && (x <= 1000) && ((x % 163) == 0));
//		if ((x >= 500) && (x <= 1000) && ((x % 163) == 0))
//			++inRange;
//	}
//
//	cout << "good numbers: " << inRange << "\n";
//
//	return 0;
//}

//int main()
//{
//	int x;
//	cout << "enter a positive number\n";
//	cin >> x;
//
//	int sumDigits = 0;
//	while (x)
//	{
//		sumDigits += (x % 10);
//		x /= 10;
//	}
//
//	cout << "sum digits: " << sumDigits << "\n";
//
//	return 0;
//}

//int main()
//{
//	int x;
//	cout << "enter a positive number\n";
//	cin >> x;
//
//	bool isPrime = true;
//	for (int i = 2; (i < x); ++i)
//		isPrime &= ((x % i) != 0);
//
//	cout << x << " is " << (isPrime? "prime" : " not prime") << "\n";
//
//	return 0;
//}

//int main()
//{
//	int x;
//	cout << "enter a positive number\n";
//	cin >> x;
//
//	bool isPrime = true;
//	for (int i = 2; (i < x) && isPrime; ++i)
//		isPrime = ((x % i) != 0);
//
//	cout << x << " is " << (isPrime ? "prime" : " not prime") << "\n";
//
//	return 0;
//}

//int main()
//{
//	int x;
//	cout << "enter a positive number\n";
//	cin >> x;
//
//	bool isPrime = true;
//	for (int i = 2; i < x; ++i)
//	{
//		if ((x % i) == 0) {
//			isPrime = false;
//			break;
//		}
//	}
//
//	cout << x << " is " << (isPrime ? "prime" : "not prime") << "\n";
//
//	return 0;
//}

//int main()
//{
//	for (int i = 1; i <= 100; ++i)
//	{
//		if (i % 13 != 0)
//			cout << i << endl;
//	}
//	return 0;
//}

//int main()
//{
//	int x;
//	cout << "enter 0-6\n";
//	cin >> x;
//
//	switch (x) {
//		case 0:
//			cout << "sunday";
//			break;
//		case 1:
//			cout << "monday";
//			break;
//		case 2: cout << "tuesday"; break;
//		case 3: cout << "wednedsday"; break;
//		case 4: cout << "thursday"; break;
//		case 5: cout << "friday"; break;
//		case 6: cout << "saturday"; break;
//		default: cout << "NOT A VALID ENTRY"; break;
//	};
//	cout << endl;
//	return 0;
//}

//int main()
//{
//	int x[10];
//	cout << x[5] << "\n";
//	cout << x[100] << "\n";
//	return 0;
//}

//int main()
//{
//	const int n = 1000;
//	double x[n];
//	for (int i = 0; i < n; ++i)
//		x[i] = static_cast<double>(rand()) / RAND_MAX;
//
//	double mi = x[0], ma = x[0], sum = 0;
//	for (int i = 0; i < n; ++i) {
//		if (x[i] < mi)
//			mi = x[i];
//		if (x[i] > ma)
//			ma = x[i];
//		sum += x[i];
//	}
//
//	cout << "min: " << mi << endl
//		<< "max: " << ma << endl
//		<< "avg: " << sum/n << endl;
//	return 0;
//}

//int main()
//{
//	unsigned primes[80000]; // over-dimensioned
//
//	unsigned nPrimesFound = 0;
//	unsigned upperBound = 100;
//
//	for (unsigned i = 2; i < upperBound; ++i) {
//		unsigned largest = static_cast<unsigned>(sqrt(i));
//		bool isPrime = true;
//		for (unsigned j = 0; j < nPrimesFound; ++j) {
//			if (primes[j] > largest)
//				break;
//			if (i % primes[j] == 0) {
//				isPrime = false;
//				break;
//			}
//		}
//		if (isPrime)
//			primes[nPrimesFound++] = i;
//	}
//
//	for (unsigned i = 0; i < nPrimesFound; ++i)
//		cout << primes[i] << ", ";
//	cout << endl;
//	return 0;
//}


int main()
{
	const int n = 100;
	int x[n];
	bool isprime = 1;
	for (int j = 0; j < n; j++)
	{
		x[j] = j + 1;
		bool isprime = 1;
		if (x[j] < 2)
			cout << x[j] <<endl;
		else
		{
			for (int i = 2; i < x[j]; i++)
			{
				if (x[j] % i == 0)
				{
					isprime = 0;

					break;
				}

			}
			if (isprime != 0)
			{
				cout << x[j] <<endl;
			}
		}
	}
}