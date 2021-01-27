// TABULAR METHOD
// The pdf table is generated offline and passed here as a buffer.
// Constansts below define the bins similar to the distribution at:
// Efficient Rendering of Heterogeneous Polydisperse Granular Media
// of Mueller et at.

#include "../Tools/Parameters.h"

// HG factor [-1,1] linear
#define BINS_G 200

// Scattering albedo [0, 0.999] linear in log(1 / (1 - alpha))
#define BINS_SA 1000

// Radius [1, 256] linear in log(r)
#define BINS_R 9

// Theta [0,pi] linear in cos(theta)
#define BINS_THETA 45

#define STRIDE_G (BINS_SA * BINS_R * BINS_THETA)
#define STRIDE_SA (BINS_R * BINS_THETA)
#define STRIDE_R (BINS_THETA)


// This includes vertex info and materials in space1, as well as acummulation buffers for a path-tracer
#include "..\Tools\CommonPT.h" 

// Random using is HybridTaus
#include "..\Tools\Randoms.h"

// Includes some functions for surface scattering and texture mapping
#include "..\Tools\Scattering.h"

#include "../Tools/HGPhaseFunction.h"

// Top level structure with the scene
RaytracingAccelerationStructure Scene : register(t0, space0);

// Shell Transport Function used as PDF of the outgoing position
// The STF is storing actually the empirical cdf(x | g, phi, r)
StructuredBuffer<float> STF						: register(t1);
// [r, phi, g] -> s^1
StructuredBuffer<float> OnceTimeScatteringAlbedo : register(t2);
//Texture3D<float> OnceTimeScatteringAlbedo		: register(t2);
// [r, phi, g] -> s^m
StructuredBuffer<float> MultipleTimeScatteringAlbedo : register(t3);
//Texture3D<float> MultipleTimeScatteringAlbedo	: register(t3);

float DistanceToSphereBoundary(float3 x, float3 w)
{
	//float a = dot(w,w); <- 1 because w is normalized
	float b = 2 * dot(x, w);
	float c = dot(x, x) - 1;

	float Disc = b * b - 4 * c;

	if (Disc <= 0)
		return 0;

	// Assuming x is inside the sphere, only the positive root is needed (intersection forward w).
	return max(0, (-b + sqrt(Disc)) / 2);
}

// Uses path tracing to generate variable outcomes.
// this allow to check the overall algorithm correctness.
bool ExactSampleCosXAndW(float g, float phi, float r, out float theta, out float beta, out float alpha)
{
	//float logR = max(0, round(log2(r)));
	//int rBin = (int)logR;
	//rBin += (random() < r / (1 << rBin) - 1); // better than interpolate beteen logR and logR+1
	//r = pow(2.0, rBin);

	float3 x = float3(0, 0, 0);
	float3 w = float3(0, 0, 1);

	theta = 1;
	beta = 0;
	alpha = 0;

	float t = -log(1 - random()) / r;
	if (t >= 1) // leaves without scattering
	{
		return true;
	}
	x += w * t; // first flight

	while (true) {
		if (random() < 1 - phi) // Absorption
			return false;

		w = ImportanceSamplePhase(g, w); // scattering event...

		float d = DistanceToSphereBoundary(x, w);

		t = -log(max(0.000000001, 1 - random())) / r;

		if (t >= d)
		{
			x += d * w;

			// code path variables in a compact way
			float3 zAxis = float3(0, 0, 1);
			float3 xAxis = abs(x.z) > 0.999 ? float3(1, 0, 0) : normalize(cross(x, float3(0, 0, 1)));
			float3 yAxis = cross(zAxis, xAxis);

			float3x3 normR = transpose(float3x3(xAxis, yAxis, zAxis));
			float3 normx = mul(x, normR);
			float3 normw = mul(w, normR);
			float3 B = float3(1, 0, 0);
			float3 T = cross(normx, B);
			theta = normx.z;
			beta = dot(normw, T);
			alpha = dot(normw, B);
			//alpha = 0;
			//beta = 0;
			//if (alpha != 0) // => N>1
			//{ // check the cosine distribution assumption...
			//	float angle = random() * 2 * 3.141596;
			//	float rad = sqrt(random());
			//	alpha = rad * sin(angle);
			//	beta = rad * cos(angle);
			//}

			return true;
		}

		x += t * w;
	}
}

bool SampleCosXAndW2(float g, float phi, float r, out float theta, out float beta, out float alpha) {
	float logR = max(0, round(log2(r)));

	int rBin = (int)logR;
	rBin += (random() < r / (1 << rBin) - 1); // better than interpolate beteen logR and logR+1
	//int rBin = (int)logR + (random() < (logR % 1)); // better than interpolate beteen logR and logR+1

	int phiBin = phi > 0.999 ? BINS_SA - 1 :
		(log(1 / (1 - phi))) * (BINS_SA - 1) / (log(1 / 0.001));
	int gBin = min((g * 0.5 + 0.5) * BINS_G, BINS_G - 1);
	int offsetInSTFTable = rBin * STRIDE_R + phiBin * STRIDE_SA + gBin * STRIDE_G;

	float selectingCase = random();
	float prob1Scat = OnceTimeScatteringAlbedo[gBin * (BINS_SA * BINS_R) + phiBin * BINS_R + rBin];
	//float prob1Scat = OnceTimeScatteringAlbedo[uint3(rBin, phiBin, gBin)];
	float probmScat = MultipleTimeScatteringAlbedo[gBin * (BINS_SA * BINS_R) + phiBin * BINS_R + rBin];
	//float probmScat = MultipleTimeScatteringAlbedo[uint3(rBin, phiBin, gBin)];

	if (selectingCase * (1 - exp(-r) - prob1Scat) < probmScat) {
		// Mueller assumes separability of p(x,w) = p(x)p(w)
		// alpha and beta are choosen uniformly in plane disc,
		// therefore they represent a cosine weighted outgoing direction

		selectingCase = random() * probmScat;

		float angle = random() * 2 * 3.141596;
		float rad = sqrt(random());
		alpha = rad * sin(angle);
		beta = rad * cos(angle);

		// theta should be sampled from the table!

		int minTheta = offsetInSTFTable;
		int maxTheta = offsetInSTFTable + BINS_THETA - 1;

		// Binary search after working
		while (minTheta < maxTheta) {
			int med = (minTheta + maxTheta) / 2;
			if (selectingCase < STF[med])
				maxTheta = med;
			else
				minTheta = med + 1;
		}
		theta = 2 * ((minTheta - offsetInSTFTable) + random()) / BINS_THETA - 1;
		return true;
	}

	alpha = 0;
	beta = 0;
	theta = 1;

	float t = -log(max(0.000000001, 1 - random())) / r;

	if (t >= 1)
		return true;

	float3 x = float3(0, 0, t);
	float3 w = ImportanceSamplePhase(g, float3(0, 0, 1));

	if (random() < 1 - phi)
		return false;

	float d = DistanceToSphereBoundary(x, w);

	t = -log(max(0.000000001, 1 - random())) / r;
	if (t < d)
		return false; // not a single scattering (=> absorption)

	x += w * d;
	float cosW = dot(x, w);
	beta = sqrt(max(0, 1 - cosW * cosW));
	theta = x.z;
	return true;
}

// Samples an outgoing position and direction using the proposal in Muller
// Return true if the sample exit, false if ray is absorbed.
bool SampleCosXAndW(float g, float phi, float r, out float theta, out float beta, out float alpha)
{
	float logR = max(0, (log2(r)));

	int rBin = (int)logR;
	//rBin += (random() < r / (1 << rBin) - 1); // better than interpolate beteen logR and logR+1
	rBin += (random() < (logR % 1)); // better than interpolate beteen logR and logR+1
	r = pow(2.0, rBin);

	int phiBin = phi > 0.999 ? BINS_SA - 1 :
		(log(1 / (1 - phi))) * (BINS_SA - 2) / (log(1 / 0.001));
	int gBin = min((g * 0.5 + 0.5) * BINS_G, BINS_G - 1);
	int offsetInSTFTable = rBin * STRIDE_R + phiBin * STRIDE_SA + gBin * STRIDE_G;

	float selectingCase = random();
	float prob0Scat = exp(-r);
	float prob1Scat = OnceTimeScatteringAlbedo[gBin*(BINS_SA * BINS_R) + phiBin * BINS_R + rBin];
	//float prob1Scat = OnceTimeScatteringAlbedo[uint3(rBin, phiBin, gBin)];
	float probmScat = MultipleTimeScatteringAlbedo[gBin * (BINS_SA * BINS_R) + phiBin * BINS_R + rBin];
	//float probmScat = MultipleTimeScatteringAlbedo[uint3(rBin, phiBin, gBin)];

	if (selectingCase < prob0Scat) // no scattering
	{
		theta = 1; // exit in straight direction
		alpha = 0;
		beta = 0; // outgoing direction is in line with the normal.
		return true;
	}
	selectingCase -= prob0Scat;
	if (selectingCase < prob1Scat) {
		// free path inside the sphere...
		// This code is Mueller proposal. It is bias, because assumes 
		// the one-time-scattering is an event that decide a free path inside and then 
		// choose independently the outgoing direction using the phase function.
		// This independece assumption is an error because the second flight is conditioned to
		// lay outside the sphere... a flight is not equally probable depending on the first decided position.
		float t = -log(1 - random() * (1 - prob0Scat)) / r;
		float3 w = ImportanceSamplePhase(g, float3(0, 0, 1)); // scattering event...
		float d = DistanceToSphereBoundary(float3(0, 0, t), w);
		float3 x = float3(0, 0, t) + w * d;
		float cosBeta = dot(x, w);
		beta = sqrt(1 - cosBeta * cosBeta);
		alpha = 0;
		theta = x.z;
		return true;

		//float cosW = w.z;
		//alpha = 0; // with only one scattering the outgoing direction is aligned to the 0.0.1 axis.
		//float sinW = sqrt(1 - cosW * cosW);
		//beta = t * sinW; // beta is the projection in plane, so really sin B
		//float sinB = beta;
		//float cosB = sqrt(1 - sinB * sinB);
		//// cos Theta = cos (W - B) = cosW cosB + sinW sinB
		//theta = cosW * cosB + sinW * sinB;
		return true;
	}
	selectingCase -= prob1Scat;
	if (selectingCase < probmScat) {
		// Mueller assumes separability of p(x,w) = p(x)p(w)
		// alpha and beta are choosen uniformly in plane disc,
		// therefore they represent a cosine weighted outgoing direction

		float angle = random() * 2 * 3.141596;
		float rad = sqrt(random());
		alpha = rad * sin(angle);
		beta = rad * cos(angle);

		// theta should be sampled from the table!

		int minTheta = offsetInSTFTable;
		int maxTheta = offsetInSTFTable + BINS_THETA - 1;

		// Binary search after working
		while (minTheta < maxTheta) {
			int med = (minTheta + maxTheta) / 2;
			if (selectingCase < STF[med])
				maxTheta = med;
			else
				minTheta = med + 1;
		}
		theta = 2 * ((minTheta - offsetInSTFTable) + random()) / BINS_THETA - 1;
		return true;

		// sequential search after working
		/*for (int thetaBin = 0; thetaBin < BINS_THETA; thetaBin++)
			if (selectingCase < STF[offsetInSTFTable + thetaBin])
			{
				theta = 2 * (thetaBin + random()) / BINS_THETA - 1;
				return true;
			}*/
	}
	alpha = beta = theta = 0;
	return false; // Absorbed ray
}

struct GridInfo {
	// Index of the base geometry (grid).
	int GridIndex;
	// Transform from the world space to the grid.
	// Considers Instance Transform, Geometry Transform and Grid transform.
	float4x4 FromWorldToGrid;
	// Scaling factor to convert distances from grid (cell is unit) to world.
	float FromGridToWorldScaling;
};
StructuredBuffer<GridInfo> GridInfos : register(t4);
Texture3D<float> DistanceField[100] : register(t5);

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

bool GenerateVariablesWithTable(float G, float Phi, float3 win, float density, out float3 x, out float3 w)
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

	float theta, beta, alpha;
	if (!SampleCosXAndW(G, Phi, density, theta, beta, alpha))
	//if (!SampleCosXAndW2(G, Phi, density, theta, beta, alpha))
	//if (!ExactSampleCosXAndW(G, Phi, density, theta, beta, alpha))
		return false;

	x = float3(0, sqrt(1 - theta * theta), theta);
	float3 N = x;
	float3 B = float3(1, 0, 0);
	float3 T = cross(x, B);
	w = normalize(N * sqrt(max(0, 1 - beta * beta - alpha * alpha)) + T * beta + B * alpha);

	x = mul(x, R);
	w = mul(w, R); // move to radial space

	return true;
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

float3 ComputePath(float3 O, float3 D, inout int complexity)
{
	int cmp = NumberOfPasses % 3;
	float3 importance = 0;
	importance[cmp] = 3;
	float3 x = O;
	float3 w = D;

	bool inMedium = false;

	int bounces = 0;

	float3 result = 0;

	bool isOutside = true;

	[loop]
	while (true)
	{
		complexity++;

		RayPayload payload = (RayPayload)0;
		if (!Intersect(x, w, payload)) // 
			return importance * (SampleSkybox(w) + SampleLight(w) * (bounces > 0));
		
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
			bool UsePT = DispatchRaysIndex().x < DispatchRaysDimensions().x* PathtracingRatio;

			if (UsePT) {
				x += t * w; // free traverse in a medium

				if (random() < 1 - volMaterial.ScatteringAlbedo[cmp]) // absorption instead
					return 0;

				w = ImportanceSamplePhase(volMaterial.G[cmp], w); // scattering event...
			}
			else
			{
				float r = MaximalRadius(x, payload.TransformIndex);

				float er = volMaterial.Extinction[cmp] * r;

				if (er >= 1)
				{
					while (er >= 1) {
						er = min(er, 256);
						float3 _x, _w;
						if (!GenerateVariablesWithTable(volMaterial.G[cmp], volMaterial.ScatteringAlbedo[cmp], w, er, _x, _w))
							return 0;
						w = _w;
						x += _x * er / volMaterial.Extinction[cmp];

						r = MaximalRadius(x, payload.TransformIndex);
						er = volMaterial.Extinction[cmp] * r;
					}
				}
				else
				{
					x += t * w; // free traverse in a medium

					if (random() < 1 - volMaterial.ScatteringAlbedo[cmp]) // absorption instead
						return 0;

					w = ImportanceSamplePhase(volMaterial.G[cmp], w); // scattering event...
				}
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

	float3 color = ComputePath(O, D, complexity);

	if (any(isnan(color)))
		color = float3(0, 0, 0);

	AccumulateOutput(raysIndex, color, complexity);
}