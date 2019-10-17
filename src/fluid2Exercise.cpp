#include "scene.h"
#include "pcg_solver.h"
#include <iostream>
#include <fstream>
#include <string>

#define EMITTER_WIDTH 10
#define EMITTER_CENTER 50
#define EMITTER_HEIGHT 5
#define INK_AMOUNT 1
#define UP_VELOCITY 10

namespace
{
	///////////////////////////////////////////////
	// Add any custom classes or functions here! //
	///////////////////////////////////////////////

	float check(Index2 & a, Array2<float> & b)
	{
		Index2 sz = b.getSize();
		if (a.x >= 0 && a.x < sz.x && a.y >= 0 && a.y < sz.y)
			return b[a];
		else
			return 0.0f;
	}

	int clamp(int x, int low, int up)
	{
		return std::max(low, std::min(x, up));
	}

}

// advection
void Fluid2::fluidAdvection(const float dt)
{
	Index2 gridSize = grid.getSize();
	Bbox2 box = grid.getDomain();
	Vec2 dmin = box.minPosition;
	Vec2 dmax = box.maxPosition;

	// ink advection
	{
		Array2<float> ink2(ink);
		const Index2& Inksize = ink.getSize();
		for (int i=0; i<gridSize.x; i++)
		{
			for (int j=0; j<gridSize.y; j++)
			{
				Index2 point(i, j);
				Index2 pointX = Index2(i+1, j);
				Index2 pointY = Index2(i, j+1);

				float velx = (velocityX[point] + velocityX[pointX]) * 0.5f;
				float vely = (velocityY[point] + velocityY[pointY]) * 0.5f;
				Vec2 v(velx, vely);
				
				Vec2 pointBef = grid.getCellPos(point) - dt * v;

				if (pointBef.x < dmin.x)
					pointBef.x = dmin.x;
				else if (pointBef.x > dmax.x)
					pointBef.x = dmax.x;

				Vec2 pointBefIndex = grid.getCellIndex(pointBef);

				Vec2 min = Vec2 (floor(pointBefIndex.x), floor(pointBefIndex.y));
				Index2 minID(min.x, min.y);
				Vec2 max = Vec2 (ceil(pointBefIndex.x), ceil(pointBefIndex.y));
				Index2 maxID(max.x, max.y);

				int minX = clamp(minID.x, 0, Inksize.x - 1);
				int maxX = clamp(maxID.x, 0, Inksize.x - 1);
				int minY = clamp(minID.y, 0, Inksize.y - 1);
				int maxY = clamp(maxID.y, 0, Inksize.y - 1);

				Index2 p1 = Index2(minX, minY);
				Index2 p2 = Index2(maxX, minY);
				Index2 p3 = Index2(minX, maxY);
				Index2 p4 = Index2(maxX, maxY);

				Vec2 delta = pointBefIndex - min;
				float alpha = delta.x;
				float beta = delta.y;
				float a = ink2[p1] * (1.0f - alpha) + ink2[p2] * alpha;
				float b = ink2[p3] * (1.0f - alpha) + ink2[p4] * alpha;
				float c = a * (1.0f - beta) + b * beta;

				ink[point] = c;
			}
		}
	}

	// velocity advection
	if (Scene::testcase >= Scene::SMOKE)
	{
		Array2<float> ink2(ink);
		const Index2& Inksize = ink.getSize();
		Array2<float> horizontalV(velocityX.getSize());
		Array2<float> verticalV(velocityY.getSize());
		Index2 faceSizeX = grid.getSizeFacesX();
		Index2 faceSizeY = grid.getSizeFacesY();

		for (int i = 0; i < faceSizeX.x; i++)
		{
			for (int j = 0; j < faceSizeX.y; j++)
			{
				Index2 point(i, j);
				Vec2 pointD = grid.getFaceXPos(point);

				float v = check(point, velocityX);

				float fi = check(point, velocityY);
				float gi = check(Index2(point.x - 1, point.y), velocityY);
				float psi = check(Index2(point.x, point.y + 1), velocityY);
				float omega = check(Index2(point.x - 1, point.y + 1), velocityY);
				float u = (((gi + fi) * 0.5f) + ((omega + psi) * 0.5f)) * 0.5f;
				
				Vec2 vel(v, u);

				Vec2 pointBef = pointD - dt * vel;

				if (pointBef.x < dmin.x)
					pointBef.x = dmin.x;
				else if (pointBef.x > dmax.x)
					pointBef.x = dmax.x;

				if (pointBef.y < dmin.y)
					pointBef.y = dmin.y;
				else if (pointBef.y > dmax.y)
					pointBef.y = dmax.y;

				Vec2 faceIndexX = grid.getFaceIndex(pointBef, 0);

				Vec2 min = Vec2(floor(faceIndexX.x), floor(faceIndexX.y));
				Index2 minID(min.x, min.y);
				Vec2 max = Vec2(ceil(faceIndexX.x), ceil(faceIndexX.y));
				Index2 maxID(max.x, max.y);

				int minX = clamp(minID.x, 0, Inksize.x - 1);
				int maxX = clamp(maxID.x, 0, Inksize.x - 1);
				int minY = clamp(minID.y, 0, Inksize.y - 1);
				int maxY = clamp(maxID.y, 0, Inksize.y - 1);

				Index2 p1 = Index2(minX, minY);
				Index2 p2 = Index2(maxX, minY);
				Index2 p3 = Index2(minX, maxY);
				Index2 p4 = Index2(maxX, maxY);

				Vec2 delta = faceIndexX - min;
				float alpha = delta.x;
				float beta = delta.y;
				float a = velocityX[p1] * (1.0f - alpha) + velocityX[p2] * alpha;
				float b = velocityX[p3] * (1.0f - alpha) + velocityX[p4] * alpha;
				float c = a * (1.0f - beta) + b * beta;

				Index2 sz = horizontalV.getSize();
				if (point.x > 0 && point.x < sz.x - 1 && point.y >= 0 && point.y < sz.y)
				{
					horizontalV[point] = c;
				}
			}
		}

		for (int i = 0; i < faceSizeY.x; i++)
		{
			for (int j = 0; j < faceSizeY.y; j++)
			{
				Index2 point(i, j);
				Vec2 pointD = grid.getFaceYPos(point);

				float fi = check(point, velocityX);
				float gi = check(Index2(point.x + 1, point.y), velocityX);
				float psi = check(Index2(point.x, point.y - 1), velocityX);
				float omega = check(Index2(point.x + 1, point.y - 1), velocityX);
				float v = (((gi + fi) * 0.5f) + ((omega + psi) * 0.5f)) * 0.5f;
				
				float u = check(point, velocityY);
				Vec2 vel(v, u);

				Vec2 pointBef = pointD - dt * vel;

				if (pointBef.x < dmin.x)
					pointBef.x = dmin.x;
				else if (pointBef.x > dmax.x)
					pointBef.x = dmax.x;

				if (pointBef.y < dmin.y)
					pointBef.y = dmin.y;
				else if (pointBef.y > dmax.y)
					pointBef.y = dmax.y;

				Vec2 faceIndexY = grid.getFaceIndex(pointBef, 1);

				Vec2 min = Vec2(floor(faceIndexY.x), floor(faceIndexY.y));
				Index2 minID(min.x, min.y);
				Vec2 max = Vec2(ceil(faceIndexY.x), ceil(faceIndexY.y));
				Index2 maxID(max.x, max.y);

				int minX = clamp(minID.x, 0, Inksize.x - 1);
				int maxX = clamp(maxID.x, 0, Inksize.x - 1);
				int minY = clamp(minID.y, 0, Inksize.y - 1);
				int maxY = clamp(maxID.y, 0, Inksize.y - 1);

				Index2 p1 = Index2(minX, minY);
				Index2 p2 = Index2(maxX, minY);
				Index2 p3 = Index2(minX, maxY);
				Index2 p4 = Index2(maxX, maxY);

				Vec2 delta = faceIndexY - min;
				float alpha = delta.x;
				float beta = delta.y;
				float a = velocityY[p1] * (1.0f - alpha) + velocityY[p2] * alpha;
				float b = velocityY[p3] * (1.0f - alpha) + velocityY[p4] * alpha;
				float c = a * (1.0f - beta) + b * beta;

				Index2 sz = verticalV.getSize();
				if (point.x >= 0 && point.x < sz.x && point.y > 0 && point.y < sz.y - 1)
				{
					verticalV[point] = c;
				}
			}
		}
		velocityX = horizontalV;
		velocityY = verticalV;
	}
}

// emission
void Fluid2::fluidEmission()
{
	if (Scene::testcase >= Scene::SMOKE)
	{

		// emit source ink
		{
			for (int i = (EMITTER_CENTER - EMITTER_WIDTH * 0.5); i < EMITTER_CENTER + EMITTER_WIDTH * 0.5; i++)
			{
				for (int j = EMITTER_HEIGHT; j < EMITTER_HEIGHT * 2; j++)
				{
					Index2 emitter(i, j);
					ink[emitter] = INK_AMOUNT;
				}
			}
		}
		// emit source velocity
		{
			for (int i = (EMITTER_CENTER - EMITTER_WIDTH * 0.5); i < EMITTER_CENTER + EMITTER_WIDTH * 0.5; i++)
			{
				for (int j = EMITTER_HEIGHT; j < EMITTER_HEIGHT * 2; j++)
				{
					Index2 emitter(i, j);
					velocityY[emitter] = UP_VELOCITY;
				}
			}
		}
	}
}

// volume forces
void Fluid2::fluidVolumeForces(const float dt)
{
	if (Scene::testcase >= Scene::SMOKE)
	{
		// gravity
		int tam = velocityY.getSize().x * velocityY.getSize().y;

		for (int i = 0; i < tam; i++)
			velocityY[i] += dt * Scene::kGravity;
	}
}

// viscosity
void Fluid2::fluidViscosity(const float dt)
{
	if (Scene::testcase >= Scene::SMOKE)
	{
		// viscosity
		Vec2 cellDx = grid.getCellDx();

		Index2 faceSizeX = grid.getSizeFacesX();
		Index2 faceSizeY = grid.getSizeFacesY();

		Array2<float> horizontalV(velocityX.getSize());
		float density = Scene::kDensity;
		float viscosity = Scene::kViscosity;
		for (int i = 0; i < faceSizeX.x; i++)
		{
			for (int j = 0; j < faceSizeX.y; j++)
			{
				Index2 p(i, j);

				//	(i,j)
				float V_ij = check(p, velocityX);
				//	(i+1,j)
				float V_right = check(Index2(p.x + 1, p.y), velocityX);
				//	(i-1,j)			
				float V_left = check(Index2(p.x - 1, p.y), velocityX);
				//	(i,j+1)
				float V_up = check(Index2(p.x, p.y+1), velocityX);
				//	(i,j-1)
				float V_down = check(Index2(p.x, p.y-1), velocityX);

				float value = V_ij + (dt / density) * viscosity * (((V_right - 2.0f * V_ij + V_left) / (cellDx.x * cellDx.x)) + ((V_up - 2.0f * V_ij + V_down) / (cellDx.y * cellDx.y)));

				Index2 sz = horizontalV.getSize();
				if (p.x > 0 && p.x < sz.x - 1 && p.y >= 0 && p.y < sz.y)
				{
					horizontalV[p] = value;
				}
			}
		}

		Array2<float> verticalV(velocityY.getSize());
		for (int i = 0; i < faceSizeY.x; i++)
		{
			for (int j = 0; j < faceSizeY.y; j++)
			{
				Index2 p(i, j);

				//	(i,j)
				float V_ij = check(p, velocityY);
				//	(i+1,j)
				float V_right = check(Index2(p.x + 1, p.y), velocityY);
				//	(i-1,j)
				float V_left = check(Index2(p.x - 1, p.y), velocityY);
				//	(i,j+1)
				float V_up = check(Index2(p.x, p.y + 1), velocityY);
				//	(i,j-1)
				float V_down = check(Index2(p.x, p.y - 1), velocityY);

				float value = V_ij + (dt / density) * viscosity * (((V_right - 2.0f * V_ij + V_left) / (cellDx.x * cellDx.x)) + ((V_up - 2.0f * V_ij + V_down) / (cellDx.y * cellDx.y)));
				
				Index2 sz = verticalV.getSize();
				if (p.x >= 0 && p.x < sz.x && p.y > 0 && p.y < sz.y - 1)
				{
					verticalV[p] = value;
				}
			}
		}
		velocityX = horizontalV;
		velocityY = verticalV;
	}
}

// pressure
void Fluid2::fluidPressureProjection(const float dt)
{
	if (Scene::testcase >= Scene::SMOKE)
	{
		// pressure
		float density = Scene::kDensity;
		double res;
		int volt;

		Index2 gridSize = grid.getSize();
		SparseMatrix<double> A(gridSize.x * gridSize.y, 5);

		const Index2& size_v = velocityY.getSize();
		for (int i = 0; i < size_v.x; i++)
			velocityY[Index2(i, 0)] = 0.0f;

		std::vector<double> den;
		den.resize(gridSize.x * gridSize.y);
		std::vector<double> b(gridSize.x * gridSize.y);

		PCGSolver<double> solver;
		solver.set_solver_parameters(1e-2, 1000);

		Vec2 cellDx = grid.getCellDx();

		for (int i = 0; i < gridSize.x; i++)
		{
			for (int j = 0; j < gridSize.y; j++)
			{
				float elementX = 2.0f / (cellDx.x * cellDx.x);
				float elementY = 2.0f / (cellDx.y * cellDx.y);
				if (i == 0 || i == gridSize.x - 1)
					elementX = 1.0f / (cellDx.x * cellDx.x);
				if (j == 0)
					elementY = 1.0f / (cellDx.y * cellDx.y);

				int linearID = pressure.getLinearIndex(i, j);
				Index2 point_ij(i, j);
				A.set_element(linearID, linearID, elementX + elementY);

				if (j == gridSize.y - 1)
					A.set_element(linearID, linearID, elementX + elementY);
			
				Index2 P_left(i - 1, j);
				if (P_left.x < gridSize.x)
					A.set_element(linearID, pressure.getLinearIndex(P_left.x, P_left.y), (-1.0f / (cellDx.x * cellDx.x)));
				Index2 P_right(i + 1, j);
				if (P_right.x < gridSize.x)
					A.set_element(linearID, pressure.getLinearIndex(P_right.x, P_right.y), (-1.0f / (cellDx.x * cellDx.x)));
				Index2 P_down(i, j - 1);
				if (P_down.y < gridSize.y)
					A.set_element(linearID, pressure.getLinearIndex(P_down.x, P_down.y), (-1.0f / (cellDx.y * cellDx.y)));
				Index2 P_up(i, j + 1);
				if (P_up.y < gridSize.y)
					A.set_element(linearID, pressure.getLinearIndex(P_up.x, P_up.y), (-1.0f / (cellDx.y * cellDx.y)));

				float v_x = check(point_ij, velocityX);
				float v_y = check(point_ij, velocityY);
				float v_right = check(P_right, velocityX);
				float v_up = check(P_up, velocityY);

				float deltaVx = (v_right - v_x) / cellDx.x;
				float deltaVy = (v_up - v_y) / cellDx.y;
				float denValue = (-density / dt) * (deltaVx + deltaVy);

				den[linearID] = denValue;
			}
		}

		solver.solve(A, den, b, res, volt);

		for (int i = 0; i < gridSize.x; i++)
		{
			for (int j = 0; j < gridSize.y; j++)
			{
				int upd = pressure.getLinearIndex(i, j);
				pressure[Index2(i, j)] = b[upd];
			}
		}

		Array2<float> velX(velocityX.getSize());
		Array2<float> velY(velocityY.getSize());
		for (int i = 0; i < gridSize.x; i++)
		{
			for (int j = 0; j < gridSize.y; j++)
			{
				Index2 p(i, j);

				float ro = check(p, pressure);

				Index2 horizontal = Index2(p.x + 1, p.y);
				float horizontalV = check(horizontal, velocityX);
				float horizontalP = check(horizontal, pressure);
				float calcX = horizontalV - (dt / density) * (horizontalP - ro) / cellDx.x;

				Index2 szX = velX.getSize();
				if (horizontal.x >= 0 && horizontal.x < szX.x && horizontal.y >= 0 && horizontal.y < szX.y)
					velX[horizontal] = calcX;

				Index2 vertical = Index2(p.x, p.y+1);
				float verticalV = check(vertical, velocityY);
				float verticalP = check(vertical, pressure);
				float calcY = verticalV - (dt / density) * (verticalP - ro) / cellDx.y;

				Index2 szY = velY.getSize();
				if (vertical.x >= 0 && vertical.x < szY.x && vertical.y >= 0 && vertical.y < szY.y)
					velY[vertical] = calcY;
			}
		}
		velocityX = velX;
		velocityY = velY;
	}
}