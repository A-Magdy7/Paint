// project1.cpp : Defines the entry point for the application.
//

#include "stdafx.h"
#include "project1.h"
#include <windows.h>
#include <Commdlg.h>
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
using namespace std;


#define MAX_LOADSTRING 100

// Global Variables:
HINSTANCE hInst;								// current instance
TCHAR szTitle[MAX_LOADSTRING];					// The title bar text
TCHAR szWindowClass[MAX_LOADSTRING];			// the main window class name

// Forward declarations of functions included in this code module:
ATOM				MyRegisterClass(HINSTANCE hInstance);
BOOL				InitInstance(HINSTANCE, int);
LRESULT CALLBACK	WndProc(HWND, UINT, WPARAM, LPARAM);
INT_PTR CALLBACK	About(HWND, UINT, WPARAM, LPARAM);

int APIENTRY _tWinMain(_In_ HINSTANCE hInstance,
                     _In_opt_ HINSTANCE hPrevInstance,
                     _In_ LPTSTR    lpCmdLine,
                     _In_ int       nCmdShow)
{
	UNREFERENCED_PARAMETER(hPrevInstance);
	UNREFERENCED_PARAMETER(lpCmdLine);

 	// TODO: Place code here.
	MSG msg;
	HACCEL hAccelTable;

	// Initialize global strings
	LoadString(hInstance, IDS_APP_TITLE, szTitle, MAX_LOADSTRING);
	LoadString(hInstance, IDC_PROJECT1, szWindowClass, MAX_LOADSTRING);
	MyRegisterClass(hInstance);

	// Perform application initialization:
	if (!InitInstance (hInstance, nCmdShow))
	{
		return FALSE;
	}

	hAccelTable = LoadAccelerators(hInstance, MAKEINTRESOURCE(IDC_PROJECT1));

	// Main message loop:
	while (GetMessage(&msg, NULL, 0, 0))
	{
		if (!TranslateAccelerator(msg.hwnd, hAccelTable, &msg))
		{
			TranslateMessage(&msg);
			DispatchMessage(&msg);
		}
	}

	return (int) msg.wParam;
}



//
//  FUNCTION: MyRegisterClass()
//
//  PURPOSE: Registers the window class.
//
ATOM MyRegisterClass(HINSTANCE hInstance)
{
	WNDCLASSEX wcex;

	wcex.cbSize = sizeof(WNDCLASSEX);

	wcex.style			= CS_HREDRAW | CS_VREDRAW;
	wcex.lpfnWndProc	= WndProc;
	wcex.cbClsExtra		= 0;
	wcex.cbWndExtra		= 0;
	wcex.hInstance		= hInstance;
	wcex.hIcon			= LoadIcon(hInstance, MAKEINTRESOURCE(IDI_PROJECT1));
	wcex.hCursor		= LoadCursor(NULL, IDC_ARROW);
	wcex.hbrBackground	= (HBRUSH)(COLOR_WINDOW+1);
	wcex.lpszMenuName	= MAKEINTRESOURCE(IDC_PROJECT1);
	wcex.lpszClassName	= szWindowClass;
	wcex.hIconSm		= LoadIcon(wcex.hInstance, MAKEINTRESOURCE(IDI_SMALL));

	return RegisterClassEx(&wcex);
}

//
//   FUNCTION: InitInstance(HINSTANCE, int)
//
//   PURPOSE: Saves instance handle and creates main window
//
//   COMMENTS:
//
//        In this function, we save the instance handle in a global variable and
//        create and display the main program window.
//
BOOL InitInstance(HINSTANCE hInstance, int nCmdShow)
{
   HWND hWnd;

   hInst = hInstance; // Store instance handle in our global variable

   hWnd = CreateWindow(szWindowClass, szTitle, WS_OVERLAPPEDWINDOW,
      CW_USEDEFAULT, 0, CW_USEDEFAULT, 0, NULL, NULL, hInstance, NULL);

   if (!hWnd)
   {
      return FALSE;
   }

   ShowWindow(hWnd, nCmdShow);
   UpdateWindow(hWnd);

   return TRUE;
}

//
//  FUNCTION: WndProc(HWND, UINT, WPARAM, LPARAM)
//
//  PURPOSE:  Processes messages for the main window.
//
//  WM_COMMAND	- process the application menu
//  WM_PAINT	- Paint the main window
//  WM_DESTROY	- post a quit message and return
//
//

enum
{
	ID_Filling, BG_White, BG_Black, BG_Red, BG_Blue, BG_Yellow, BG_Green,
	CIRCLE_CAR, CIRCLE_POLAR, CIRCLE_IPOLAR, CIRCLE_MP,
	CURVE_FD, CURVE_SD, CURVE_HER, CURVE_BEZ, CURVE_SP,
	FILE_NEW, FILE_SAVE, FILE_OPEN, FILE_Exit, 
	LINE_DIRECT, LINE_DDA, LINE_MP,
	REC_CLIP_P, REC_CLIP_L,
	CIR_CLIP_P, CIR_CLIP_L
};

const COLORREF black = RGB(0, 0, 0), white = RGB(255, 255, 255), blue = RGB(0, 0, 255), red = RGB(255,0,0), green = RGB(0,255,0);
COLORREF DColor = black, WColor = white;

void createMenu(HWND hwnd)
{
	HMENU hMenubar = CreateMenu();
	HMENU hFile = CreateMenu();
	HMENU hBackground = CreateMenu();
	HMENU hLine = CreateMenu();
	HMENU hCircle = CreateMenu();
	HMENU hCurve = CreateMenu();
	HMENU hClipping = CreateMenu();
	HMENU hThirdDegree = CreateMenu();
	HMENU hCRectangle = CreateMenu();
	HMENU hCCircle = CreateMenu();

	AppendMenu(hMenubar, MF_POPUP, (UINT_PTR)hFile, _T("File"));
	AppendMenu(hMenubar, MF_POPUP, (UINT_PTR)hLine, _T("Line"));
	AppendMenu(hMenubar, MF_POPUP, (UINT_PTR)hCircle, _T("Circle"));
	AppendMenu(hMenubar, MF_POPUP, (UINT_PTR)hCurve, _T("Curve"));
	AppendMenu(hMenubar, MF_STRING, ID_Filling, _T("Convex Filling"));
	AppendMenu(hMenubar, MF_POPUP, (UINT_PTR)hClipping, _T("Clipping"));

	AppendMenu(hClipping, MF_POPUP, (UINT_PTR)hCRectangle, _T("Rectangle"));
	AppendMenu(hClipping, MF_POPUP, (UINT_PTR)hCCircle, _T("Circle"));

	AppendMenu(hCRectangle, MF_STRING, REC_CLIP_P, _T("Point Clipping"));
	AppendMenu(hCRectangle, MF_STRING, REC_CLIP_L, _T("Line Clipping"));

	AppendMenu(hCCircle, MF_STRING, CIR_CLIP_P, _T("Point Clipping"));
	AppendMenu(hCCircle, MF_STRING, CIR_CLIP_L, _T("Line Clipping"));

	AppendMenu(hFile, MF_STRING, FILE_NEW, _T("New"));
	AppendMenu(hFile, MF_STRING, FILE_SAVE, _T("Save"));
	AppendMenu(hFile, MF_STRING, FILE_OPEN, _T("Open"));
	AppendMenu(hFile, MF_POPUP, (UINT_PTR)hBackground, _T("Back Ground Color"));
	AppendMenu(hFile, MF_STRING, FILE_Exit, _T("Exit"));

	AppendMenu(hBackground, MF_STRING, BG_White, _T("white"));
	AppendMenu(hBackground, MF_STRING, BG_Black, _T("black"));
	AppendMenu(hBackground, MF_STRING, BG_Red, _T("red"));
	AppendMenu(hBackground, MF_STRING, BG_Blue, _T("blue"));
	AppendMenu(hBackground, MF_STRING, BG_Green, _T("green"));

	AppendMenu(hLine, MF_STRING, LINE_DIRECT, _T("Direct"));
	AppendMenu(hLine, MF_STRING, LINE_DDA, _T("DDA"));
	AppendMenu(hLine, MF_STRING, LINE_MP, _T("Mid Point"));

	AppendMenu(hCircle, MF_STRING, CIRCLE_CAR, _T("Cartesian"));
	AppendMenu(hCircle, MF_STRING, CIRCLE_POLAR, _T("Polar"));
	AppendMenu(hCircle, MF_STRING, CIRCLE_IPOLAR, _T("Iterative Polar"));
	AppendMenu(hCircle, MF_STRING, CIRCLE_MP, _T("Mid Point"));

	AppendMenu(hCurve, MF_STRING, CURVE_FD, _T("First Degree"));
	AppendMenu(hCurve, MF_STRING, CURVE_SD, _T("Second Degree"));
	AppendMenu(hCurve, MF_POPUP, (UINT_PTR)hThirdDegree, _T("Third Degree"));

	AppendMenu(hThirdDegree, MF_STRING, CURVE_HER, _T("Hermite"));
	AppendMenu(hThirdDegree, MF_STRING, CURVE_BEZ, _T("Bezier"));
	AppendMenu(hThirdDegree, MF_STRING, CURVE_SP, _T("Splines"));

	SetMenu(hwnd, hMenubar);
}

void ChangeBGColor(HDC hdc, HWND hwnd, COLORREF col)
{
	RECT rect;
	HBRUSH brush = CreateSolidBrush(col);

	GetWindowRect(hwnd, &rect);
	FillRect(hdc, &rect, brush);

	if (col == white)
		DColor = black;
	else
		DColor = white;
}

string Browse(bool f)
{
	char filename[MAX_PATH];

	OPENFILENAMEA ofn;
	ZeroMemory(&filename, sizeof(filename));
	ZeroMemory(&ofn, sizeof(ofn));
	ofn.lStructSize = sizeof(ofn);
	ofn.hwndOwner = NULL;  // If you have a window to center over, put its HANDLE here
	ofn.lpstrFilter = "Text Files\0*.txt\0Any File\0*.*\0";
	ofn.lpstrFile = filename;
	ofn.nMaxFile = MAX_PATH;
	ofn.lpstrTitle = "Select File";
	ofn.Flags = (f ? OFN_FILEMUSTEXIST : NULL);
	if (GetOpenFileNameA(&ofn))
		return string(filename);
	else
		return "";
}

void Save(HDC hdc, HWND hwnd)
{
	string Fname = Browse(0);

	if (Fname == "")
		return;

	RECT rect;
	COLORREF col;
	ofstream out;

	GetWindowRect(hwnd, &rect);
	out.open(Fname);

	for (int i = rect.left; i < rect.right; i++)
		for (int j = rect.top; j < rect.bottom; j++)
			col = GetPixel(hdc, i, j), out.write((char*)&col, sizeof(col));

	out.close();
}

void Open(HDC hdc, HWND hwnd)
{
	string Fname = Browse(1);

	if (Fname == "")
		return;

	RECT rect;
	COLORREF col;
	ifstream in;

	GetWindowRect(hwnd, &rect);
	in.open(Fname);

	for (int i = rect.left; i < rect.right; i++)
		for (int j = rect.top; j < rect.bottom; j++)
			in.read((char*)&col, sizeof(col)), SetPixel(hdc, i, j, col);

	in.close();
}

void swap(int &x, int &y)
{
	int tmp = x;
	x = y;
	y = tmp;
}

struct point{
	int x, y;
	point(int a, int b){
		x = a;
		y = a;
	}
	point(){
		x = 0;
		y = 0;
	}
};

void Draw8Points(HDC hdc, int xc, int yc, int a, int b, COLORREF col)
{
	SetPixel(hdc, xc + a, yc + b, col);
	SetPixel(hdc, xc + a, yc - b, col);
	SetPixel(hdc, xc - a, yc + b, col);
	SetPixel(hdc, xc - a, yc - b, col);
	SetPixel(hdc, xc + b, yc + a, col);
	SetPixel(hdc, xc + b, yc - a, col);
	SetPixel(hdc, xc - b, yc + a, col);
	SetPixel(hdc, xc - b, yc - a, col);
}

void cartesian(HDC hdc, int xc, int yc, int R, COLORREF col)
{
	int x = 0, y = R;

	Draw8Points(hdc, xc, yc, x, y, col);

	while (x <= y)
	{
		x++;
		y = round(sqrt(R*R - x*x));
		Draw8Points(hdc, xc, yc, x, y, col);
	}
}

void polar(HDC hdc, int xc, int yc, int R, COLORREF col)
{
	int x = R, y = 0;
	double s = 0, ds = 1.0 / R;

	Draw8Points(hdc, xc, yc, x, y, col);

	while (x >= y)
	{
		s += ds;
		x = round(R*cos(s));
		y = round(R*sin(s));

		Draw8Points(hdc, xc, yc, x, y, col);
	}
}

void IterativePolar(HDC hdc, int xc, int yc, int R, COLORREF col)
{
	double x = R, y = 0, cs = cos(1.0 / R), ss = sin(1.0 / R), tmp;

	Draw8Points(hdc, xc, yc, x, y, col);

	while (x >= y)
	{
		tmp = x*cs - y*ss;

		y = x*ss + y*cs;
		x = tmp;

		Draw8Points(hdc, xc, yc, round(x), round(y), col);
	}
}

void CircleMP(HDC hdc, int xc, int yc, int R, COLORREF col)
{
	int x = 0, y = R, d = 1 - R;

	Draw8Points(hdc, xc, yc, x, y, col);

	while (x < y)
	{
		if (d < 0)
			d += 2 * x + 3;
		else
			y--, d += 2 * (x - y) + 5;

		x++;

		Draw8Points(hdc, xc, yc, x, y, col);
	}
}

void DirectLine(HDC hdc, int xs, int ys, int xe, int ye, COLORREF col)
{
	int dx = xe - xs, dy = ye - ys, x, y;
	double s;

	if (abs(dx) > abs(dy))
	{
		if (xe < xs)
			swap(xe, xs), swap(ye, ys);

		s = dy*1.0 / dx, x = xs;

		SetPixel(hdc, xs, ys, col);

		while (x < xe)
		{
			y = round(ys + (x - xs)*s);

			SetPixel(hdc, ++x, y, col);
		}
	}

	else
	{
		if (ye < ys)
			swap(xe, xs), swap(ye, ys);

		s = dx*1.0 / dy, y = ys;

		SetPixel(hdc, xs, ys, col);

		while (y < ye)
		{
			x = round(xs + (y - ys)*s);

			SetPixel(hdc, x, ++y, col);
		}
	}
}

void DDA(HDC hdc, int xs, int ys, int xe, int ye, COLORREF col)
{
	int steps, dx = xe - xs, dy = ye - ys;
	double x=xs, y=ys, Xinc, Yinc;

	steps = max(abs(dx), abs(dy));

	Xinc = dx*1.0 / steps;
	Yinc = dy*1.0 / steps;

	SetPixel(hdc, x, y, col);

	while (steps--)
	{
		x += Xinc;
		y += Yinc;

		SetPixel(hdc, x, y, col);
	}
}

void LineMP(HDC hdc, int xs, int ys, int xe, int ye, COLORREF col)
{
	int dx = abs(xe-xs), dy = abs(ye-ys), d, c1, c2, Xinc, Yinc, x, y;

	if (dx>dy)
	{
		if (xe < xs)
			swap(xe, xs), swap(ye, ys);

		Yinc = ye > ys ? 1 : -1;
		d = dx - 2 * dy;
		c1 = 2 * (dx - dy);
		c2 = -2 * dy;
		x = xs, y = ys;

		SetPixel(hdc, x, y, col);

		while (x < xe)
		{
			if (d < 0)
				y += Yinc, d += c1;
			else
				d += c2;

			SetPixel(hdc, ++x, y, col);
		}
	}

	else
	{
		if (ye < ys)
			swap(xe, xs), swap(ye, ys);

		Xinc = xe > xs ? 1 : -1;
		d = 2 * dx - dy;
		c1 = 2 * (dx - dy);
		c2 = 2 * dx;
		x = xs, y = ys;

		SetPixel(hdc, x, y, col);

		while (y < ye)
		{
			if (d > 0)
				x += Xinc, d += c1;
			else
				d += c2;

			SetPixel(hdc, x, ++y, col);
		}
	}
}

/********************************************* First Degree *****************************************************/
void first_degree(HDC hdc, point ps, point pe, COLORREF color){
	typedef int Beta, Alpha;
	Beta B1 = ps.x;
	Beta B2 = ps.y;
	Alpha A1 = pe.x - ps.x;
	Alpha A2 = pe.y - ps.y;
	int n = max(abs(A1), abs(A2));
	//int n = 100;
	double dt = 1.0 / n;
	double x = (double)ps.x;
	double y = (double)ps.y;
	double t = 0;
	for (int i = 0; i < n; i++){
		SetPixel(hdc, round(x), round(y), color);
		x = A1 *t + B1;
		y = A2 *t + B2;
		t += dt;
	}
}

void second_degree(HDC hdc, point p1, point T1, point p2, COLORREF color){
	typedef int Alpha;
	typedef int Beta;
	typedef int Gama;
	Beta B1, B2;
	B1 = 2*(T1.x-p1.x);
	B2 = 2*(T1.y-p1.y);
	Alpha A1 = p2.x + p1.x - 2*T1.x;
	Alpha A2 = p2.y + p1.y - 2*T1.y;
	Gama G1, G2;
	G1 = p1.x;
	G2 = p1.y;
	int n = max(abs(p2.x - p1.x), abs(p2.y - p1.y));
	double dt = 1.0 / (n - 1);
	int x = p1.x;
	int y = p1.y;
	double t = 0;
	for (int i = 0; i < n; i++){
		SetPixel(hdc, round(x), round(y), color);
		double t2 = t*t;
		x = A1 * t2 + B1 * t + G1;
		y = A2 * t2 + B2 * t + G2;
		t += dt;
	}
}

typedef double vec[4];
typedef double matrix[4][4];
void mul(matrix A, vec B, vec C){
	for (int i = 0; i < 4; i++){
		C[i] = 0;
		for (int j = 0; j < 4; j++){
			C[i] += A[i][j] * B[j];
		}
	}
}

double dot(vec a, vec b){
	double sum = 0;
	for (int i = 0; i < 4; i++){
		sum += a[i] * b[i];
	}
	return sum;
}

void DrawCurveHermit(HDC hdc, point p1, point t1, point p2, point t2, COLORREF color){
	static matrix H = { { 2, 1, -2, 1 }, { -3, -2, 3, -1 }, { 0, 1, 0, 0 }, { 1, 0, 0, 0 } };
	vec Vx = { p1.x, t1.x, p2.x, t2.x };
	vec Vy = { p1.y, t1.y, p2.y, t2.y };
	vec Gx, Gy;
	mul(H, Vx, Gx);
	mul(H, Vy, Gy);
	//int n = max(abs(p2.x - p1.x), abs(p2.y - p1.y));
	int n = 10000;
	double dt = 1.0 / n;
	double t = 0;
	for (int i = 0; i < n; i++){
		double t2 = t*t;
		double t3 = t2*t;
		vec Vt = { t3, t2, t, 1 };
		double x = dot(Gx, Vt);
		double y = dot(Gy, Vt);
		SetPixel(hdc, round(x), round(y), color);
		t += dt;
	}
}
/****************************************************************************************************************/

/********************************************** Bezeir Curve ****************************************************/

void DrawCurveBezeir(HDC hdc, point p1, point p2, point p3, point p4, COLORREF color)
{
	point T1(3 * (p2.x - p1.x), 3 * (p2.y - p1.y));
	point T2(3 * (p4.x - p3.x), 3 * (p4.y - p3.y));
	DrawCurveHermit(hdc, p1, T1, p2, T2, color);
}

/****************************************************************************************************************/
void DrawCardinalSpline(HDC hdc, point P[], int n, double c, COLORREF color)
{
	double c1 = c / 2;
	point T0(c1*(P[2].x - P[0].x), c1*(P[2].y - P[0].y));
	second_degree(hdc, P[0], T0, P[1], color);
	for (int i = 2; i<n - 1; i++)
	{
		point T1(c1*(P[i + 1].x - P[i - 1].x), c1*(P[i + 1].y - P[i - 1].y));
		DrawCurveHermit(hdc, P[i - 1], T0, P[i], T1, color);
		T0 = T1;
	}
	second_degree(hdc, P[n - 1], T0, P[n - 2], color);
}

struct entry {
	int xleft; // xmin
	int xright;  // xmax
	entry() {
		xleft = 800;
		xright = 0;
	}
	entry(int a, int b) {
		xleft = a;
		xright = b;
	}
};

void initial_entries(entry table[], int n) {  // utility function used to initialize the entry table
	for (int i = 0; i < n; i++) {
		table[i].xleft = INT_MAX;  // 800 is considered as the max value
		table[i].xright = INT_MIN;   // 0 is considered as the min value
	}
}

void ScanEdge(point v1, point v2, entry table[]) {
	//takes the end points of the edge, generates all points of intersection with scan lines and update the entry table
	if (v1.y == v2.y)return;
	if (v2.y < v1.y) {
		// swap

		int temp = v1.x;
		v1.x = v2.x;
		v2.x = temp;

		temp = v1.y;
		v1.y = v2.y;
		v2.y = temp;
	}
	double mi = double(v2.x - v1.x) / (v2.y - v1.y);  // inverse of slope
	int y = v1.y;
	double x = (double)v1.x;
	while (y < v2.y) {
		if (x < table[y].xleft)
			table[y].xleft = (int)ceil(x); // xmin
		if (x > table[y].xright)
			table[y].xright = (int)floor(x);  // xmax
		y++;
		x += mi;
	}
}

void DrawLines(HDC hdc, entry tbl[], COLORREF color) {
	// utility to draw scan lines stored in the entry table.
	for (int i = 0; i < 1000; i++)
	{
		if ((tbl[i].xleft < tbl[i].xright))
			DirectLine(hdc, tbl[i].xleft, i, tbl[i].xright, i, color);
	}
}

void convexFill(HDC hdc, point p[]/*vertices points*/, int n, COLORREF color) {
	// the main function of the algorithm that contains the overall logic of the algorithm
	entry * table = new entry[1000];
	initial_entries(table, 1000);
	point v1 = p[n - 1];
	for (int i = 0; i < n; i++) {
		point v2 = p[i];
		ScanEdge(v1, v2, table);
		v1 = p[i];
	}
	DrawLines(hdc, table, color);
	delete []table;
}


/*********************************************************************************************/
void DrawPolygon(HDC hdc, point p[], int n, COLORREF color) {
	point v1 = p[n - 1];
	for (int i = 0; i < n; i++) {
		point v2 = p[i];
		DirectLine(hdc, v1.x, v1.y, v2.x, v2.y, color);
		v1 = p[i];
	}
}

/**********************************************************************************************/
void PointClippingtocircle(int x, int y, int xs, int ys, double R, HDC hdc, COLORREF color)
{
	double r = sqrt((double)(x - xs)*(x - xs) + (y - ys)*(y - ys));
	if (r <= R) SetPixel(hdc, x, y, color);
}

bool check(int x, int y, int xs, int ys, double R)
{
	double r = sqrt((double)(x - xs)*(x - xs) + (y - ys)*(y - ys));
	if (r <= R) return true;
	else return false;
}

void drawLine(HDC hdc, int xs, int ys, int xe, int ye, int xc, int yc, double R, COLORREF color) {

	int dx = xe - xs;
	int dy = ye - ys;
	SetPixel(hdc, xs, ys, color);
	if (abs(dx) >= abs(dy))
	{
		int x = xs, xinc = dx>0 ? 1 : -1;
		double y = ys, yinc = (double)dy / dx*xinc;
		while (x != xe)
		{
			x += xinc;
			y += yinc;
			if (check(x, y, xc, yc, R)){

				SetPixel(hdc, x, round(y), color);
			}
		}
	}
	else
	{
		int y = ys, yinc = dy>0 ? 1 : -1;
		double x = xs, xinc = (double)dx / dy*yinc;
		while (y != ye)
		{
			x += xinc;
			y += yinc;
			if (check(x, y, xc, yc, R)){

				SetPixel(hdc, round(x), y, color);
			}
		}
	}
}
/********************************************************************************************************/
void PointClippingtoRec(int x, int y, int Xmin, int Xmax, int Ymin, int Ymax, HDC hdc, COLORREF color)
{
	if (x>=Xmin && x<=Xmax && y<=Ymax && y>=Ymin)
		SetPixel(hdc, x, y, color);
}

void DrawRec(HDC hdc, int xs, int xe, int ys, int ye, COLORREF color)
{
	LineMP(hdc, xs, ys, xs, ye, color);
	LineMP(hdc, xe, ys, xe, ye, color);

	LineMP(hdc, xe, ye, xs, ye, color);
	LineMP(hdc, xe, ys, xs, ys, color);
}

union Code
{
	unsigned All : 4;
	struct{ unsigned left : 1, top : 1, right : 1, bottom : 1; };
};

Code getCode(double x, double y, int Xmin, int Xmax, int Ymin, int Ymax)
{
	Code code;

	code.All = 0;

	if (x < Xmin)
		code.left = 1;
	else if (x > Xmax)
		code.right = 1;
	if (y < Ymin)
		code.top = 1;
	else if (y > Ymax)
		code.bottom = 1;

	return code;
}

void YInter(double xs, double ys, double xe, double ye, int x, double &xi, double &yi)
{
	xi = x;
	yi = ys + (x - xs)*(ye - ys) / (xe - xs);
}

void XInter(double xs, double ys, double xe, double ye, int y, double &xi, double &yi)
{
	yi = y;
	xi = xs + (y - ys)*(xe - xs) / (ye - ys);
}

void LineClippingtoRec(double xs, double ys, double xe, double ye, int Xmin, int Xmax, int Ymin, int Ymax, HDC hdc, COLORREF color)
{
	Code code1 = getCode(xs, ys, Xmin, Xmax, Ymin, Ymax);
	Code code2 = getCode(xe, ye, Xmin, Xmax, Ymin, Ymax);

	while ((code1.All || code2.All) && !(code1.All & code2.All))
	{
		double xi, yi;

		if (code1.All)
		{
			if (code1.left)
				YInter(xs, ys, xe, ye, Xmin, xi, yi);

			else if (code1.right)
				YInter(xs, ys, xe, ye, Xmax, xi, yi);

			else if (code1.top)
				XInter(xs, ys, xe, ye, Ymin, xi, yi);

			else
				XInter(xs, ys, xe, ye, Ymax, xi, yi);

			xs = xi, ys = yi;

			code1 = getCode(xs, ys, Xmin, Xmax, Ymin, Ymax);
		}

		else if (code2.All)
		{
			if (code2.left)
				YInter(xs, ys, xe, ye, Xmin, xi, yi);

			else if (code2.right)
				YInter(xs, ys, xe, ye, Xmax, xi, yi);

			else if (code2.top)
				XInter(xs, ys, xe, ye, Ymin, xi, yi);

			else
				XInter(xs, ys, xe, ye, Ymax, xi, yi);

			xe = xi, ye = yi;

			code2 = getCode(xe, ye, Xmin, Xmax, Ymin, Ymax);
		}
	}

	if (!code1.All && !code2.All)
		LineMP(hdc, xs, ys, xe, ye, color);
}

void init(HWND hWnd)
{
	RECT rect;
	GetWindowRect(hWnd, &rect);
	InvalidateRect(hWnd, &rect, false);
}

LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
	int wmId, wmEvent;
	PAINTSTRUCT ps;
	HDC hdc = GetDC(hWnd);
	static int funID, count=7, sz;
	static point p[5];

	switch (message)
	{
	case WM_CREATE:
		createMenu(hWnd);
		break;
	case WM_COMMAND:
		wmId    = LOWORD(wParam);
		wmEvent = HIWORD(wParam);
		// Parse the menu selections:
		switch (wmId)
		{
		case FILE_SAVE:
			Save(hdc,hWnd);
			break;
		case FILE_OPEN:
			Open(hdc,hWnd);
			break;
		case FILE_Exit:
			DestroyWindow(hWnd);
		case FILE_NEW:
		case BG_White:
			count = -1, funID = 0, WColor = white, init(hWnd);
			break;
		case BG_Black:
			count = -1, funID = 0, WColor = black, init(hWnd);
			break;
		case BG_Red:
			count = -1, funID = 0, WColor = red, init(hWnd);
			break;
		case BG_Green:
			count = -1, funID = 0, WColor = green, init(hWnd);
			break;
		case BG_Blue:
			count = -1, funID = 0, WColor = blue, init(hWnd);
			break;
		case ID_Filling:
			funID = 1, count=4, sz=4;
			break;
		case LINE_DIRECT:
			funID = 2, count = 1, sz = 1;
			break;
		case LINE_DDA:
			funID = 3, count = 1, sz = 1;
			break;
		case LINE_MP:
			funID = 4, count = 1, sz = 1;
			break;
		case CIRCLE_CAR:
			funID = 5, count = 1, sz = 1;
			break;
		case CIRCLE_POLAR:
			funID = 6, count = 1, sz = 1;
			break;
		case CIRCLE_IPOLAR:
			funID = 7, count = 1, sz = 1;
			break;
		case CIRCLE_MP:
			funID = 8, count = 1, sz = 1;
			break;
		case CURVE_FD:
			funID = 9, count = 1, sz = 1;
			break;
		case CURVE_SD:
			funID = 10, count = 2, sz = 2;
			break;
		case CURVE_HER:
			funID = 11, count = 3, sz = 3;
			break;
		case CURVE_BEZ:
			funID = 12, count = 3, sz = 3;
			break;
		case CURVE_SP:
			funID = 13, count = 4, sz = 4;
			break;
		case CIR_CLIP_P:
			funID = 14, count = 2, sz = 2;
			break;
		case CIR_CLIP_L:
			funID = 15, count = 3, sz = 3;
			break;
		case REC_CLIP_P:
			funID = 16, count = 2, sz = 2;
			break;
		case REC_CLIP_L:
			funID = 17, count = 3, sz = 3;
			break;
		case IDM_ABOUT:
			DialogBox(hInst, MAKEINTRESOURCE(IDD_ABOUTBOX), hWnd, About);
			break;
		case IDM_EXIT:
			DestroyWindow(hWnd);
			break;
		default:
			return DefWindowProc(hWnd, message, wParam, lParam);
		}
		break;
	case WM_LBUTTONDOWN:
		if (count == 0)
		{
			p[sz - count].x = LOWORD(lParam);
			p[sz - count].y = HIWORD(lParam);
			count--;
			init(hWnd);
		}
		else if (count == 1)
		{
			p[sz - count].x = LOWORD(lParam);
			p[sz - count].y = HIWORD(lParam);
			count--;
		}
		else if (count == 2)
		{
			p[sz - count].x = LOWORD(lParam);
			p[sz - count].y = HIWORD(lParam);
			count--;
		}
		else if (count == 3)
		{
			p[sz - count].x = LOWORD(lParam);
			p[sz - count].y = HIWORD(lParam);
			count--;
		}
		else if (count == 4)
		{
			p[sz - count].x = LOWORD(lParam);
			p[sz - count].y = HIWORD(lParam);
			count--;
		}
		break;
	case WM_PAINT:
		if (count == -1)
		{
			hdc = BeginPaint(hWnd, &ps);

			if (funID == 0)
			{
				ChangeBGColor(hdc, hWnd, WColor);
				count = 7;
			}

			if (funID == 1)
			{
				DrawPolygon(hdc, p, 5, DColor);
				convexFill(hdc, p, 5, DColor);
				count = 7;
			}

			else if (funID == 2)
			{
				DirectLine(hdc, p[0].x, p[0].y, p[1].x, p[1].y, DColor);
				count = 7;
			}

			else if (funID == 3)
			{
				DDA(hdc, p[0].x, p[0].y, p[1].x, p[1].y, DColor);
				count = 7;
			}

			else if (funID == 4)
			{
				LineMP(hdc, p[0].x, p[0].y, p[1].x, p[1].y, DColor);
				count = 7;
			}

			else if (funID == 5)
			{
				double r = sqrt((p[1].x - p[0].x)*(p[1].x - p[0].x) + (p[1].y - p[0].y)*(p[1].y - p[0].y));
				cartesian(hdc, p[0].x, p[0].y, r, DColor);
				count = 7;
			}

			else if (funID == 6)
			{
				double r = sqrt((p[1].x - p[0].x)*(p[1].x - p[0].x) + (p[1].y - p[0].y)*(p[1].y - p[0].y));
				polar(hdc, p[0].x, p[0].y, r, DColor);
				count = 7;
			}

			else if (funID == 7)
			{
				double r = sqrt((p[1].x - p[0].x)*(p[1].x - p[0].x) + (p[1].y - p[0].y)*(p[1].y - p[0].y));
				IterativePolar(hdc, p[0].x, p[0].y, r, DColor);
				count = 7;
			}

			else if (funID == 8)
			{
				double r = sqrt((p[1].x - p[0].x)*(p[1].x - p[0].x) + (p[1].y - p[0].y)*(p[1].y - p[0].y));
				CircleMP(hdc, p[0].x, p[0].y, r, DColor);
				count = 7;
			}

			else if (funID == 9)
			{
				first_degree(hdc, p[0], p[1], DColor);
				count = 7;
			}

			else if (funID == 10)
			{
				second_degree(hdc, p[0], p[1], p[2], DColor);
				count = 7;
			}

			else if (funID == 11)
			{
				DrawCurveHermit(hdc, p[0], p[1], p[2], p[3], DColor);
				count = 7;
			}

			else if (funID == 12)
			{
				DrawCurveBezeir(hdc, p[0], p[1], p[2], p[3], DColor);
				count = 7;
			}

			else if (funID == 13)
			{
				DrawCardinalSpline(hdc, p, sz, 5, DColor);
				count = 7;
			}

			else if (funID == 14)
			{
				double r = sqrt((p[1].x - p[0].x)*(p[1].x - p[0].x) + (p[1].y - p[0].y)*(p[1].y - p[0].y));
				CircleMP(hdc, p[0].x, p[0].y, r, DColor);
				PointClippingtocircle(p[2].x, p[2].y, p[0].x, p[0].y, r, hdc, DColor);
				count = 7;
			}

			else if (funID == 15)
			{
				double r = sqrt((p[1].x - p[0].x)*(p[1].x - p[0].x) + (p[1].y - p[0].y)*(p[1].y - p[0].y));
				CircleMP(hdc, p[0].x, p[0].y, r, DColor);
				drawLine(hdc, p[2].x, p[2].y, p[3].x, p[3].y, p[0].x, p[0].y, r, DColor);
				count = 7;
			}

			else if (funID == 16)
			{
				if (p[0].x > p[1].x)
					swap(p[0].x, p[1].x);

				if (p[0].y > p[1].y)
					swap(p[0].y, p[1].y);

				DrawRec(hdc, p[0].x, p[1].x, p[0].y, p[1].y, DColor);

				PointClippingtoRec(p[2].x, p[2].y, p[0].x, p[1].x, p[0].y, p[1].y, hdc, DColor);
				count = 7;
			}

			else if (funID == 17)
			{
				if (p[0].x > p[1].x)
					swap(p[0].x, p[1].x);

				if (p[0].y > p[1].y)
					swap(p[0].y, p[1].y);

				DrawRec(hdc, p[0].x, p[1].x, p[0].y, p[1].y, DColor);

				LineClippingtoRec(p[2].x, p[2].y, p[3].x, p[3].y, p[0].x, p[1].x, p[0].y, p[1].y, hdc, DColor);
				count = 7;
			}

			EndPaint(hWnd, &ps);
		}
		break;
	case WM_DESTROY:
		PostQuitMessage(0);
		break;
	default:
		return DefWindowProc(hWnd, message, wParam, lParam);
	}
	return 0;
}

// Message handler for about box.
INT_PTR CALLBACK About(HWND hDlg, UINT message, WPARAM wParam, LPARAM lParam)
{
	UNREFERENCED_PARAMETER(lParam);
	switch (message)
	{
	case WM_INITDIALOG:
		return (INT_PTR)TRUE;

	case WM_COMMAND:
		if (LOWORD(wParam) == IDOK || LOWORD(wParam) == IDCANCEL)
		{
			EndDialog(hDlg, LOWORD(wParam));
			return (INT_PTR)TRUE;
		}
		break;
	}
	return (INT_PTR)FALSE;
}
