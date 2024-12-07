#include <math.h>
#include <stdio.h>
#define Y1(x) (x * x -3*x+4)
#define Y2(x) (x * x * x * x + 5 * x * x *x+ 3*x*x-4)
#define Y3(x) (sin(2*x) + cos(x*x))
int main()
{
double start_value = 0, end_value = 3, tolerance = 0.4;
double y1[30], y2[30], y3[30];
int count = 0;
for (double temp = start_value; temp <= end_value; temp += tolerance) {
y1[count] = Y1(temp);
y2[count] = Y2(temp);
y3[count] = Y3(temp);
count++;
}
printf("\nX\n");
for (double temp = start_value; temp <= end_value; temp += tolerance) {
printf("%.4lf | ", temp);
}
printf("\n\nF(1)\n");
for (int i = 0; i < count; i++) {
printf("%.4lf | ", y1[i]);
}
printf("\n\nF(2)\n");
for (int i = 0; i < count; i++) {
printf("%.4lf | ", y2[i]);
}
printf("\n\nF(3)\n");
for (int i = 0; i < count; i++) {
printf("%.4lf | ", y3[i]);

}
return 0;
}
