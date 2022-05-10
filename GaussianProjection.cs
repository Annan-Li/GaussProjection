using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using OpenCvSharp;

namespace draw
{
	internal class GaussianProjection
	{
		static double FE = 500000;
		static double FN = 0;
		static double k0 = 1.0;

		static double r = 6371393;						// 地球平均半径

		static double a = 6378137;						// 长半轴 国家2000坐标系
		static double f = 1.0 / 298.257222101;			// 扁率

		//static double a = 6378245;					// 长半轴 北京54坐标系
		//static double f = 1 / 298.3;					// 扁率

		/// <summary>
		/// 高斯投影正算
		/// </summary>
		/// <param name="p">经纬度坐标，X:纬度，Y:经度，单位：度</param>
		/// <param name="CalMeridian">中央子午线经度，单位：度</param>
		/// <returns>高斯投影坐标（南北，东西），单位：m</returns>
		public static Point2d Projection(Point2d p ,double CalMeridian)
		{
			double e, ep, B, L, B0, L0, k1, k2, k3, k4, T, C, A, N, M, M0, X, Y;

			e= Math.Sqrt(2.0 * (double)f - (double)f * (double)f);
			ep = Math.Sqrt((1.0 / (double)((1.0 - f) * (1.0 - f))) - 1.0);                      // 参考椭球第二偏心率
			B = (double)((double)(p.X) / 180.0) * Math.PI;										// 计算点纬度（弧度）
			L = (double)((double)(p.Y) / 180.0) * Math.PI;										// 计算点经度（弧度）
			B0 = 0;                                                                             // 原点纬度（弧度）
			L0 = CalMeridian / 180.0 * Math.PI;                                                 // 中央子午线经度（弧度）

			k1 = 1.0 - e * e / 4.0 - 3.0 * e * e * e * e / 64.0 - 5.0 * e * e * e * e * e * e / 256.0;
			k2 = -3.0 * e * e / 8.0 - 3.0 * e * e * e * e / 32.0 - 45.0 * e * e * e * e * e * e / 1024.0;
			k3 = 15.0 * e * e * e * e / 256.0 + 45.0 * e * e * e * e * e * e / 1024.0;
			k4 = -35.0 * e * e * e * e * e * e / 3072.0;

			T = Math.Tan(B) * Math.Tan(B);
			C = ep * ep * Math.Cos(B) * Math.Cos(B);
			A = (L - L0) * Math.Cos(B);
			N = a / (Math.Sqrt(1 - e * e * Math.Sin(B) * Math.Sin(B)));
			M = a * (k1 * B + k2 * Math.Sin(2.0 * B) + k3 * Math.Sin(4.0 * B) + k4 * Math.Sin(6.0 * B));
			M0 = a * (k1 * B0 + k2 * Math.Sin(2.0 * B0) + k3 * Math.Sin(4.0 * B0) + k4 * Math.Sin(6.0 * B0));

			// 南北方向
			X = FN + k0 * (M - M0 + N * Math.Tan(B) * (A * A / 2.0 + (5.0 - T + 9.0 * C + 4.0 * C * C) * A * A * A * A / 24.0) + (61.0 - 58.0 * T + T * T + 600.0 * C - 330.0 * ep * ep) * A * A * A * A * A * A / 720.0);
			// 东西方向
			Y = FE + k0 * N * (A + (1 - T + C) * A * A * A / 6.0 + (5.0 - 18.0 * T + T * T + 72.0 * C - 58.0 * ep * ep) * A * A * A * A * A / 120.0);
			
			Point2d Gp = new Point2d (X, Y);
			
			return Gp;
		}

		/// <summary>
		/// 高斯投影坐标反解
		/// </summary>
		/// <param name="p">高斯投影坐标，X:南北方向，Y:东西方向，单位：m</param>
		/// <param name="CalMeridian">中央子午线经度，单位：度</param>
		/// <returns>经纬度坐标（纬度，经度），单位：度</returns>
		public static Point2d InverseProjection(Point2d p, double CalMeridian)
		{
			double L0, X, Y, e, ep, ef, Mf, B1, Bf, Nf, Rf, Tf, Cf, D, L, B;

			L0 = CalMeridian / 180.0 * Math.PI;								// 中央子午线经度（弧度）
			X = p.X;
			Y = p.Y;

			e = Math.Sqrt(2.0 * (double)f - (double)f * (double)f);
			ep = Math.Sqrt((1.0 / (double)((1.0 - f) * (1.0 - f))) - 1.0);

			ef = f / (2.0 - f);

			Mf = (X - FN) / k0;

			B1 = Mf / (a * (1.0 - e * e / 4.0 - 3.0 * e * e * e * e / 64.0 - 5.0 * e * e * e * e * e * e / 256.0));

			Bf = B1 + (3.0 * ef / 2.0 - 27.0 * ef * ef * ef / 32.0) * Math.Sin(2.0 * B1) + (21.0 * ef * ef / 16.0 - 55.0 * ef * ef * ef * ef / 32.0) * Math.Sin(4.0 * B1) + 151.0 * ef * ef * ef / 96.0 * Math.Sin(6.0 * B1) + 1097.0 * ef * ef * ef * ef / 512.0 * Math.Sin(8.0 * B1);

			Nf = a / Math.Sqrt(1 - e * e * Math.Sin(Bf) * Math.Sin(Bf));
			Rf = a * (1 - e * e) / Math.Pow((1 - e * e * Math.Sin(Bf) * Math.Sin(Bf)), 1.5);
			Tf = Math.Tan(Bf) * Math.Tan(Bf);
			Cf = ep * ep * Math.Cos(Bf) * Math.Cos(Bf);
			D = (Y - FE) / (Nf * k0);

			// 经度，单位：度
			L = (L0 + (1 / Math.Cos(Bf)) * (D - (1.0 + 2.0 * Tf + Cf) * D * D * D / 6.0 + (5.0 - 2.0 * Cf + 28.0 * Tf - 3.0 * Cf * Cf + 8.0 * ep * ep + 24.0 * Tf * Tf) * D * D * D * D * D / 120.0)) / Math.PI * 180;
			// 纬度，单位：度
			B = (Bf - Nf * Math.Tan(Bf) / Rf * (D * D / 2.0 - (5.0 + 3.0 * Tf + 10.0 * Cf - 4.0 * Cf * Cf - 9.0 * ep * ep) * D * D * D * D / 24.0 + (61.0 + 90.0 * Tf + 298.0 * Cf + 45.0 * Tf * Tf - 252.0 * ep * ep - 3.0 * Cf * Cf) * D * D * D * D * D * D / 720.0)) / Math.PI * 180;

			Point2d llp = new Point2d(B, L);

			return llp;
		}
		/// <summary>
		/// 两高斯投影坐标点计算距离，并且进行了高程修正
		/// </summary>
		/// <param name="p1">高斯投影坐标点1，单位：m</param>
		/// <param name="p2">高斯投影坐标点2，单位：m</param>
		/// <param name="H">高程，单位：m</param>
		/// <returns>两点间距离，单位：m</returns>
		public static double ProjectionDistance(Point2d p1,Point2d p2,double H)
		{
			double Distance;
			Distance = Math.Sqrt(Math.Pow(p1.X - p2.X, 2) + Math.Pow(p1.Y - p2.Y, 2));

			// 修正高程带来的距离误差
			Distance = Distance * (r + H) / r;

			return Distance;
		}
	}
}
