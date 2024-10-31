// perftesthelpers.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>

extern void maintest_rsqrt(void);
extern void maintest_sndmix(void);
extern void maintest_diffusecolor(void);
extern void maintest_dotproduct(void);
extern void maintest_lerpmeshvertexes(void);
extern void maintest_projectdlighttexture(void);

int main()
{
    printf("Hello Test!\n");

    printf("\nRun this from command line if you want the perf results!\nDebug target is bad for SSE, and Local Windows Debbuger messes with Release.\n\n");

#ifdef _DEBUG
    printf("*** Debug Mode! ***\n\n");
#endif

    //maintest_rsqrt();
    //maintest_sndmix();
    maintest_dotproduct();

    printf("\n");
    system("pause");
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
