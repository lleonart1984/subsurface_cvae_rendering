#include "../Tools/Parameters.h"

// This includes vertex info and materials in space1, as well as acummulation buffers for a path-tracer
#include "..\Tools\CommonPT.h" 

// Random using is HybridTaus
#include "..\Tools\Randoms.h"

// Includes some functions for surface scattering and texture mapping
#include "..\Tools\Scattering.h"

// Top level structure with the scene
RaytracingAccelerationStructure Scene : register(t0, space0);

struct GridInfo {
	// Index of the base geometry (grid).
	int GridIndex;
	// Transform from the world space to the grid.
	// Considers Instance Transform, Geometry Transform and Grid transform.
	float4x4 FromWorldToGrid;
	// Scaling factor to convert distances from grid (cell is unit) to world.
	float FromGridToWorldScaling;
};
StructuredBuffer<GridInfo> GridInfos : register(t1);
Texture3D<float> DistanceField[100] : register(t2);

/// Query the distance field grid.
float MaximalRadius(float3 P, int object) {

	GridInfo info = GridInfos[object];
	float3 positionInGrid = mul(float4(P, 1), info.FromWorldToGrid).xyz;
	float radius = DistanceField[info.GridIndex][positionInGrid];

	if (radius < 0) // no empty cell
		return 0;

	float3 distToMinCorner = positionInGrid % 1;
	float3 m = min(distToMinCorner, 1 - distToMinCorner);
	float minDistanceToCellBorder = min(m.x, min(m.y, m.z));
	//float safeDistanceInGridSpace = radius;
	float safeDistanceInGridSpace = minDistanceToCellBorder + radius;
	return safeDistanceInGridSpace * info.FromGridToWorldScaling;
}

cbuffer Lighting : register(b0) {
	float3 LightPosition;
	float3 LightIntensity;
	float3 LightDirection;
}

#include "../Tools/CommonEnvironment.h"

// CB used to transform rays from projection space to world.
cbuffer Transforming : register(b1) {
	float4x4 FromProjectionToWorld;
}

#include "../Tools/HGPhaseFunction.h"

//#include "CVAEScatteringModel.h"
#include "CVAEScatteringModelX.h"

float sampleNormal(float mu, float logVar) {
	//return mu + gauss() * exp(logVar * 0.5);
	return mu + gauss() * exp(clamp(logVar, -16, 16) * 0.5);
}

bool GenerateVariablesWithModel(float G, float Phi, float3 win, float density, out float3 x, out float3 w)
{
	x = float3(0, 0, 0);
	w = win;

	float3 temp = abs(win.x) >= 0.9999 ? float3(0, 0, 1) : float3(1, 0, 0);
	float3 winY = normalize(cross(temp, win));
	float3 winX = cross(win, winY);
	float rAlpha = random() * 2 * pi;
	float3x3 R = (mul(float3x3(
		cos(rAlpha), -sin(rAlpha), 0,
		sin(rAlpha), cos(rAlpha), 0,
		0, 0, 1), float3x3(winX, winY, win)));

	float codedDensity = density;// pow(density / 400.0, 0.125);

	float2 lenLatent = randomStdNormal2();
	// Generate length
	float lenInput[4];
	float lenOutput[2];
	lenInput[0] = codedDensity;
	lenInput[1] = G;
	lenInput[2] = lenLatent.x;
	lenInput[3] = lenLatent.y;
	lenModel(lenInput, lenOutput);

	float logN = max(0, sampleNormal(lenOutput[0], lenOutput[1]));
	float n = (exp(logN));
	//logN = log(n);

	if (random() >= pow(Phi, n))
		return false;

	float4 pathLatent14 = randomStdNormal4();
	float pathLatent5 = randomStdNormal();
	// Generate path
	float pathInput[8];
	float pathOutput[6];
	pathInput[0] = codedDensity;
	pathInput[1] = G;
	pathInput[2] = logN;
	pathInput[3] = pathLatent14.x;
	pathInput[4] = pathLatent14.y;
	pathInput[5] = pathLatent14.z;
	pathInput[6] = pathLatent14.w;
	pathInput[7] = pathLatent5.x;
	pathModel(pathInput, pathOutput);
	float3 sampling = randomStdNormal3();
	float3 pathMu = float3(pathOutput[0], pathOutput[1], pathOutput[2]);
	float3 pathLogVar = float3(pathOutput[3], pathOutput[4], pathOutput[5]);
	float3 pathOut = clamp(pathMu + exp(clamp(pathLogVar, -16, 16) * 0.5) * sampling, -0.9999, 0.9999);
	float costheta = pathOut.x;
	float wt = n > 1 ? pathOut.y : 0.0; // only if n > 1
	float wb = n > 2 ? pathOut.z : 0.0; // only if n > 2
	x = float3(0, sqrt(1 - costheta * costheta), costheta);
	float3 N = x;
	float3 B = float3(1, 0, 0);
	float3 T = cross(x, B);

	w = normalize(N * sqrt(max(0, 1 - wt * wt - wb * wb)) + T * wt + B * wb);
	x = mul(x, (R));
	w = mul(w, (R)); // move to radial space

	return true;// random() >= 1 - pow(Phi, n);
}

bool GenerateFullVariablesWithModel(float G, float Phi, float3 win, float density, out float3 x, out float3 w, out float3 X, out float3 W, out float factor)
{
	x = float3(0, 0, 0);
	w = win;
	X = x;
	W = w;
	factor = 1;

	float3 temp = abs(win.x) >= 0.9999 ? float3(0, 0, 1) : float3(1, 0, 0);
	float3 winY = normalize(cross(temp, win));
	float3 winX = cross(win, winY);
	float rAlpha = random() * 2 * pi;
	float3x3 R = (mul(float3x3(
		cos(rAlpha), -sin(rAlpha), 0,
		sin(rAlpha), cos(rAlpha), 0,
		0, 0, 1), float3x3(winX, winY, win)));

	float codedDensity = density;// pow(density / 400.0, 0.125);

	float2 lenLatent = randomStdNormal2();
	// Generate length
	float lenInput[4];
	float lenOutput[2];
	lenInput[0] = codedDensity;
	lenInput[1] = G;
	lenInput[2] = lenLatent.x;
	lenInput[3] = lenLatent.y;
	lenModel(lenInput, lenOutput);

	float logN = max(0, sampleNormal(lenOutput[0], lenOutput[1]));
	float n = round(exp(logN)+0.49);
	logN = log(n);

	if (random() >= pow(Phi, n))
		return false;

	float4 pathLatent14 = randomStdNormal4();
	float pathLatent5 = randomStdNormal();
	// Generate path
	float pathInput[8];
	float pathOutput[6];
	pathInput[0] = codedDensity;
	pathInput[1] = G;
	pathInput[2] = logN;
	pathInput[3] = pathLatent14.x;
	pathInput[4] = pathLatent14.y;
	pathInput[5] = pathLatent14.z;
	pathInput[6] = pathLatent14.w;
	pathInput[7] = pathLatent5.x;
	pathModel(pathInput, pathOutput);
	float3 sampling = randomStdNormal3();
	float3 pathMu = float3(pathOutput[0], pathOutput[1], pathOutput[2]);
	float3 pathLogVar = float3(pathOutput[3], pathOutput[4], pathOutput[5]);
	float3 pathOut = clamp(pathMu + exp(clamp(pathLogVar, -16, 16) * 0.5) * sampling, -0.9999, 0.9999);
	float costheta = pathOut.x;
	float wt = n >= 2 ? pathOut.y : 0.0; // only if n > 1
	float wb = n >= 3 ? pathOut.z : 0.0; // only if n > 2
	x = float3(0, sqrt(1 - costheta * costheta), costheta);
	float3 N = x;
	float3 B = float3(1, 0, 0);
	float3 T = cross(x, B);

	w = normalize(N * sqrt(max(0, 1 - wt * wt - wb * wb)) + T * wt + B * wb);
	x = mul(x, (R));
	w = mul(w, (R)); // move to radial space

	// Generate Scattering
	float4 scatLatent14 = randomStdNormal4();
	float scatLatent5 = randomStdNormal();
	float scatInput[12];
	float scatOutput[12];
	scatInput[0] = codedDensity;
	scatInput[1] = G;
	scatInput[2] = pow(1 - Phi, 1.0/6.0);
	scatInput[3] = logN;
	scatInput[4] = costheta;
	scatInput[5] = wt;
	scatInput[6] = wb;
	scatInput[7] = scatLatent14.x;
	scatInput[8] = scatLatent14.y;
	scatInput[9] = scatLatent14.z;
	scatInput[10] = scatLatent14.w;
	scatInput[11] = scatLatent5;

	scatModel(scatInput, scatOutput);

	if (n >= 2) {
		X = float3(
			sampleNormal(scatOutput[0], scatOutput[6]),
			sampleNormal(scatOutput[1], scatOutput[7]),
			sampleNormal(scatOutput[2], scatOutput[8]));
		X /= max(1, length(X));
		W = normalize(float3(
			sampleNormal(scatOutput[3], scatOutput[9]),
			sampleNormal(scatOutput[4], scatOutput[10]),
			sampleNormal(scatOutput[5], scatOutput[11])));
		//W /= max(0.000001, length(W));
		X = mul(X, (R));
		W = mul(W, (R)); // move to radial space

		float accum = Phi >= 0.99999 ? n : Phi * (1 - pow(Phi, n)) / (1 - Phi);
		factor = accum / pow(Phi, n);
	}

	return true;// random() >= 1 - pow(Phi, n);
}


struct RayPayload // Only used for raycasting
{
	int TriangleIndex;
	int VertexOffset;
	int MaterialIndex;
	int TransformIndex;
	float3 Barycentric;
};

/// Use RTX TraceRay to detect a single intersection. No recursion is necessary
/// Path-tracing is implemented iteratively, no recursively.
bool Intersect(float3 P, float3 D, out RayPayload payload) {
	payload = (RayPayload)0;
	RayDesc ray;
	ray.Origin = P;
	ray.Direction = D;
	ray.TMin = 0.0;
	ray.TMax = 100.0;
	TraceRay(Scene, RAY_FLAG_FORCE_OPAQUE, 0xFF, 0, 1, 0, ray, payload);
	return payload.TriangleIndex >= 0;
}

// Represents a single bounce of path tracing
// Will accumulate emissive and direct lighting modulated by the carrying importance
// Will update importance with scattered ratio divided by pdf
// Will output scattered ray to continue with
void SurfelScattering(inout float3 x, inout float3 w, inout float3 importance, Vertex surfel, Material material)
{
	float3 V = -w;

	float NdotV;
	bool invertNormal;
	float3 fN;
	float4 R, T;
	ComputeImpulses(V, surfel, material,
		NdotV,
		invertNormal,
		fN,
		R,
		T);

	float3 ratio;
	float3 direction;
	float pdf;
	RandomScatterRay(V, fN, R, T, material, ratio, direction, pdf);

	// Update gathered Importance to the viewer
	importance *= max(0, ratio);
	// Update scattered ray
	w = direction;
	x = surfel.P + sign(dot(direction, fN)) * 0.001 * fN;
}

float3 ComputePathWithNEE(float3 O, float3 D, inout int complexity)
{
	int cmp = NumberOfPasses % 3;
	float3 importance = 0;
	importance[cmp] = 3;
	float3 x = O;
	float3 w = D;

	bool inMedium = false;

	int bounces = 0;

	float3 directContribution = 0;

	bool isOutside = true;

	[loop]
	while (true)
	{
		complexity++;

		RayPayload payload = (RayPayload)0;
		if (!Intersect(x, w, payload)) // 
			return directContribution + importance * (SampleSkybox(w) + SampleLight(w) * (complexity <= 2));

		Vertex surfel = (Vertex)0;
		Material material = (Material)0;
		VolumeMaterial volMaterial = (VolumeMaterial)0;
		GetHitInfo(
			payload.Barycentric,
			payload.MaterialIndex,
			payload.TriangleIndex,
			payload.VertexOffset,
			payload.TransformIndex,
			surfel, material, volMaterial, 0, 0);

		float d = length(surfel.P - x); // Distance to the hit position.
		float t = isOutside || volMaterial.Extinction[cmp] == 0 ? 100000000 : -log(max(0.000000000001, 1 - random())) / volMaterial.Extinction[cmp];

		[branch]
		if (t >= d)
		{
			bounces += isOutside;

			if (bounces >= MAX_PATHTRACING_BOUNCES)
				return 0;

			SurfelScattering(x, w, importance, surfel, material);

			if (any(material.Specular) && material.Roulette.w > 0)
				isOutside = dot(surfel.N, w) >= 0;
		}
		else
		{ // Volume scattering or absorption
			x += t * w; // free traverse in a medium

			bool UsePT = DispatchRaysIndex().x < DispatchRaysDimensions().x* PathtracingRatio;

			if (UsePT)
			{
				if (random() < 1 - volMaterial.ScatteringAlbedo[cmp]) // absorption instead
					return directContribution;

				// Ray cast to light source
				// NEE using straight line as a naive approximation missing refraction at interface
				payload = (RayPayload)0;
				if (!Intersect(x, LightDirection, payload)) // 
					return directContribution; // BAD GEOMETRY!!
				GetHitInfo(
					payload.Barycentric,
					payload.MaterialIndex,
					payload.TriangleIndex,
					payload.VertexOffset,
					payload.TransformIndex,
					surfel, material, volMaterial, 0, 0);

				float d = length(surfel.P - x); // Get medium traversing distance

				payload = (RayPayload)0;
				bool clearPathToLight = !Intersect(surfel.P + LightDirection * 0.01, LightDirection, payload); // Shadow ray

				directContribution += clearPathToLight * importance * LightIntensity * exp(-d * volMaterial.Extinction[cmp]) * EvalPhase(volMaterial.G[cmp], w, LightDirection) / 2;

				w = ImportanceSamplePhase(volMaterial.G[cmp], w); // scattering event...
			}
			else
			{
				//return GetColor(payload.TransformIndex);
				float r = MaximalRadius(x, payload.TransformIndex);
				float er = volMaterial.Extinction[cmp] * r;

				//if (er >= 1)
				//{
				float3 _x, _w, _X, _W;
				float factor;
				if (!GenerateFullVariablesWithModel(volMaterial.G[cmp], volMaterial.ScatteringAlbedo[cmp], w, er, _x, _w, _X, _W, factor))
					return directContribution; // case of absorption

				// Gets the sample scattering position for direct lighting
				float3 xs = x + _X * r;

				// Ray cast to light source
				// NEE using straight line as a naive approximation missing refraction at interface
				payload = (RayPayload)0;
				if (!Intersect(xs, LightDirection, payload)) // 
					return directContribution; // BAD GEOMETRY!!
				GetHitInfo(
					payload.Barycentric,
					payload.MaterialIndex,
					payload.TriangleIndex,
					payload.VertexOffset,
					payload.TransformIndex,
					surfel, material, volMaterial, 0, 0);

				float d = length(surfel.P - xs); // Get medium traversing distance

				payload = (RayPayload)0;
				bool clearPathToLight = !Intersect(surfel.P + LightDirection * 0.01, LightDirection, payload); // Shadow ray

				directContribution += factor * clearPathToLight * importance * LightIntensity * exp(-d * volMaterial.Extinction[cmp]) * EvalPhase(volMaterial.G[cmp], _W, LightDirection) / 2;

				w = _w;
				x += _x * r;
				/*}
				else { // Use single scattering
					if (random() < 1 - volMaterial.ScatteringAlbedo[cmp]) // absorption instead
						return 0;
					w = ImportanceSamplePhase(volMaterial.G[cmp], w); // scattering event...
				}*/
			}
		}
	}
}

[shader("closesthit")]
void OnClosestHit(inout RayPayload payload, in BuiltInTriangleIntersectionAttributes attr)
{
	payload.Barycentric = float3(1 - attr.barycentrics.x - attr.barycentrics.y, attr.barycentrics.x, attr.barycentrics.y);
	GetIndices(payload.TransformIndex, payload.MaterialIndex, payload.TriangleIndex, payload.VertexOffset);
}

[shader("miss")]
void OnMiss(inout RayPayload payload)
{
	payload.TriangleIndex = -1;
}

[shader("raygeneration")]
void RayGen()
{
	uint2 raysIndex = DispatchRaysIndex().xy;

	uint2 raysDimensions = DispatchRaysDimensions().xy;
	StartRandomSeedForRay(raysDimensions, 1, raysIndex, 0, NumberOfPasses);

	float2 coord = (raysIndex.xy + float2(random(), random())) / raysDimensions;

	float4 ndcP = float4(2 * coord - 1, 0, 1);
	ndcP.y *= -1;
	float4 ndcT = ndcP + float4(0, 0, 1, 0);

	float4 viewP = mul(ndcP, FromProjectionToWorld);
	viewP.xyz /= viewP.w;
	float4 viewT = mul(ndcT, FromProjectionToWorld);
	viewT.xyz /= viewT.w;

	float3 O = viewP.xyz;
	float3 D = normalize(viewT.xyz - viewP.xyz);

	int complexity = 0;

	float3 color = ComputePathWithNEE(O, D, complexity);

	if (any(isnan(color)))
		color = float3(0, 0, 0);

	AccumulateOutput(raysIndex, color, complexity);
}