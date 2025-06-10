#include "Matrix.h"
#include "MathUtils.h"

Matrix4x4 Multiply(const Matrix4x4& m1, const Matrix4x4& m2) {
	Matrix4x4 result = {};

	for (int row = 0; row < 4; ++row) {
		for (int col = 0; col < 4; ++col) {
			result.m[row][col] = 0.0f;
			for (int k = 0; k < 4; ++k) {
				result.m[row][col] += m1.m[row][k] * m2.m[k][col];
			}
		}
	}

	return result;
}

Matrix4x4 Matrix4x4::operator*(const Matrix4x4& other) const
{
	return Multiply(*this, other);
}

Matrix4x4 MakeIdentityMatrix() {
	Matrix4x4 result{};
	result.m[0][0] = 1.0f;
	result.m[1][1] = 1.0f;
	result.m[2][2] = 1.0f;
	result.m[3][3] = 1.0f;
	return result;
}

// 逆行列
Matrix4x4 Inverse(const Matrix4x4& m)
{
	Matrix4x4 result{};
	float inv[16], det;
	const float* src = &m.m[0][0];

	inv[0] = src[5] * src[10] * src[15] -
		src[5] * src[11] * src[14] -
		src[9] * src[6] * src[15] +
		src[9] * src[7] * src[14] +
		src[13] * src[6] * src[11] -
		src[13] * src[7] * src[10];

	inv[4] = -src[4] * src[10] * src[15] +
		src[4] * src[11] * src[14] +
		src[8] * src[6] * src[15] -
		src[8] * src[7] * src[14] -
		src[12] * src[6] * src[11] +
		src[12] * src[7] * src[10];

	inv[8] = src[4] * src[9] * src[15] -
		src[4] * src[11] * src[13] -
		src[8] * src[5] * src[15] +
		src[8] * src[7] * src[13] +
		src[12] * src[5] * src[11] -
		src[12] * src[7] * src[9];

	inv[12] = -src[4] * src[9] * src[14] +
		src[4] * src[10] * src[13] +
		src[8] * src[5] * src[14] -
		src[8] * src[6] * src[13] -
		src[12] * src[5] * src[10] +
		src[12] * src[6] * src[9];

	inv[1] = -src[1] * src[10] * src[15] +
		src[1] * src[11] * src[14] +
		src[9] * src[2] * src[15] -
		src[9] * src[3] * src[14] -
		src[13] * src[2] * src[11] +
		src[13] * src[3] * src[10];

	inv[5] = src[0] * src[10] * src[15] -
		src[0] * src[11] * src[14] -
		src[8] * src[2] * src[15] +
		src[8] * src[3] * src[14] +
		src[12] * src[2] * src[11] -
		src[12] * src[3] * src[10];

	inv[9] = -src[0] * src[9] * src[15] +
		src[0] * src[11] * src[13] +
		src[8] * src[1] * src[15] -
		src[8] * src[3] * src[13] -
		src[12] * src[1] * src[11] +
		src[12] * src[3] * src[9];

	inv[13] = src[0] * src[9] * src[14] -
		src[0] * src[10] * src[13] -
		src[8] * src[1] * src[14] +
		src[8] * src[2] * src[13] +
		src[12] * src[1] * src[10] -
		src[12] * src[2] * src[9];

	inv[2] = src[1] * src[6] * src[15] -
		src[1] * src[7] * src[14] -
		src[5] * src[2] * src[15] +
		src[5] * src[3] * src[14] +
		src[13] * src[2] * src[7] -
		src[13] * src[3] * src[6];

	inv[6] = -src[0] * src[6] * src[15] +
		src[0] * src[7] * src[14] +
		src[4] * src[2] * src[15] -
		src[4] * src[3] * src[14] -
		src[12] * src[2] * src[7] +
		src[12] * src[3] * src[6];

	inv[10] = src[0] * src[5] * src[15] -
		src[0] * src[7] * src[13] -
		src[4] * src[1] * src[15] +
		src[4] * src[3] * src[13] +
		src[12] * src[1] * src[7] -
		src[12] * src[3] * src[5];

	inv[14] = -src[0] * src[5] * src[14] +
		src[0] * src[6] * src[13] +
		src[4] * src[1] * src[14] -
		src[4] * src[2] * src[13] -
		src[12] * src[1] * src[6] +
		src[12] * src[2] * src[5];

	inv[3] = -src[1] * src[6] * src[11] +
		src[1] * src[7] * src[10] +
		src[5] * src[2] * src[11] -
		src[5] * src[3] * src[10] -
		src[9] * src[2] * src[7] +
		src[9] * src[3] * src[6];

	inv[7] = src[0] * src[6] * src[11] -
		src[0] * src[7] * src[10] -
		src[4] * src[2] * src[11] +
		src[4] * src[3] * src[10] +
		src[8] * src[2] * src[7] -
		src[8] * src[3] * src[6];

	inv[11] = -src[0] * src[5] * src[11] +
		src[0] * src[7] * src[9] +
		src[4] * src[1] * src[11] -
		src[4] * src[3] * src[9] -
		src[8] * src[1] * src[7] +
		src[8] * src[3] * src[5];

	inv[15] = src[0] * src[5] * src[10] -
		src[0] * src[6] * src[9] -
		src[4] * src[1] * src[10] +
		src[4] * src[2] * src[9] +
		src[8] * src[1] * src[6] -
		src[8] * src[2] * src[5];

	det = src[0] * inv[0] + src[1] * inv[4] + src[2] * inv[8] + src[3] * inv[12];

	det = 1.0f / det;

	Matrix4x4 resultMatrix;
	for (int i = 0; i < 16; ++i)
	{
		((float*)resultMatrix.m)[i] = inv[i] * det;
	}

	return resultMatrix;
}


// 平行移動行列
Matrix4x4 MakeTranslateMatrix(const Vector3& translate)
{
	Matrix4x4 matrix = {};

	// 単位行列
	for (int i = 0; i < 4; ++i)
	{
		matrix.m[i][i] = 1.0f;
	}

	// 平行移動成分を設定
	matrix.m[3][0] = translate.x;
	matrix.m[3][1] = translate.y;
	matrix.m[3][2] = translate.z;

	return matrix;
}

// 拡大縮小行列
Matrix4x4 MakeScaleMatrix(const Vector3& scale)
{
	Matrix4x4 matrix = {};

	// 拡大縮小成分を設定
	matrix.m[0][0] = scale.x;
	matrix.m[1][1] = scale.y;
	matrix.m[2][2] = scale.z;
	matrix.m[3][3] = 1.0f;

	return matrix;
}

// X軸回転行列
Matrix4x4 MakeRotateXMatrix(float radian)
{
	Matrix4x4 matrix = {};

	matrix.m[0][0] = 1.0f;
	matrix.m[1][1] = std::cosf(radian);
	matrix.m[1][2] = std::sinf(radian);
	matrix.m[2][1] = -std::sinf(radian);
	matrix.m[2][2] = std::cosf(radian);
	matrix.m[3][3] = 1.0f;

	return matrix;
}

// Y軸回転行列
Matrix4x4 MakeRotateYMatrix(float radian)
{
	Matrix4x4 matrix = {};

	matrix.m[0][0] = std::cosf(radian);
	matrix.m[0][2] = -std::sinf(radian);
	matrix.m[1][1] = 1.0f;
	matrix.m[2][0] = std::sinf(radian);
	matrix.m[2][2] = std::cosf(radian);
	matrix.m[3][3] = 1.0f;

	return matrix;
}

// Z軸回転行列
Matrix4x4 MakeRotateZMatrix(float radian)
{
	Matrix4x4 matrix = {};

	matrix.m[0][0] = std::cosf(radian);
	matrix.m[0][1] = std::sinf(radian);
	matrix.m[1][0] = -std::sinf(radian);
	matrix.m[1][1] = std::cosf(radian);
	matrix.m[2][2] = 1.0f;
	matrix.m[3][3] = 1.0f;

	return matrix;
}

// 座標変換
Vector3 Transform(const Vector3& vector, const Matrix4x4& matrix)
{
	float x = vector.x * matrix.m[0][0] + vector.y * matrix.m[1][0] + vector.z * matrix.m[2][0] + matrix.m[3][0];
	float y = vector.x * matrix.m[0][1] + vector.y * matrix.m[1][1] + vector.z * matrix.m[2][1] + matrix.m[3][1];
	float z = vector.x * matrix.m[0][2] + vector.y * matrix.m[1][2] + vector.z * matrix.m[2][2] + matrix.m[3][2];
	float w = vector.x * matrix.m[0][3] + vector.y * matrix.m[1][3] + vector.z * matrix.m[2][3] + matrix.m[3][3];

	if (w != 0.0f && w != 1.0f) {
		x /= w;
		y /= w;
		z /= w;
	}

	return { x, y, z };
}

// アフィン変換行列
Matrix4x4 MakeAffineMatrix(const Vector3& scale, const Vector3& rotate, const Vector3& translate)
{
	Matrix4x4 scaleMatrix = MakeScaleMatrix(scale);
	Matrix4x4 rotateXMatrix = MakeRotateXMatrix(rotate.x);
	Matrix4x4 rotateYMatrix = MakeRotateYMatrix(rotate.y);
	Matrix4x4 rotateZMatrix = MakeRotateZMatrix(rotate.z);
	Matrix4x4 translateMatrix = MakeTranslateMatrix(translate);

	// 回転行列を合成
	Matrix4x4 rotateMatrix = Multiply(Multiply(rotateXMatrix, rotateYMatrix), rotateZMatrix);

	// 最終的なアフィン行列
	Matrix4x4 affineMatrix = Multiply(Multiply(scaleMatrix, rotateMatrix), translateMatrix);

	return affineMatrix;
}

// 正射影行列
Matrix4x4 MakeOrthographicMatrix(float left, float top, float right, float bottom, float nearClip, float farClip)
{
	Matrix4x4 matrix = {};

	matrix.m[0][0] = 2.0f / (right - left);
	matrix.m[1][1] = 2.0f / (top - bottom);
	matrix.m[2][2] = 1.0f / (farClip - nearClip);
	matrix.m[3][0] = (right + left) / (left - right);
	matrix.m[3][1] = (top + bottom) / (bottom - top);
	matrix.m[3][2] = nearClip / (nearClip - farClip);
	matrix.m[3][3] = 1.0f;

	return matrix;
}

//透視投影行列
Matrix4x4 MakePerspectiveFovMatrix(float fovY, float aspectRatio, float nearClip, float farClip)
{
	Matrix4x4 matrix = {};
	float f = 1.0f / std::tanf(fovY * 0.5f);
	matrix.m[0][0] = f / aspectRatio;
	matrix.m[1][1] = f;
	matrix.m[2][2] = farClip / (farClip - nearClip);
	matrix.m[2][3] = 1.0f;
	matrix.m[3][2] = (-farClip * nearClip) / (farClip - nearClip);
	return matrix;
}

// ビューポート変換行列
Matrix4x4 MakeViewportMatrix(float left, float top, float width, float height, float minDepth, float maxDepth)
{
	Matrix4x4 matrix = {};

	matrix.m[0][0] = width / 2.0f;
	matrix.m[1][1] = -height / 2.0f;
	matrix.m[2][2] = maxDepth - minDepth;
	matrix.m[3][0] = left + width / 2.0f;
	matrix.m[3][1] = top + height / 2.0f;
	matrix.m[3][2] = minDepth;
	matrix.m[3][3] = 1.0f;

	return matrix;
}

Vector3 WorldToScreen(const Vector3& worldPos, const Matrix4x4& WVPMatrix, float screenWidth, float screenHeight) {
	Vector3 clipPos = Transform(worldPos, WVPMatrix);

	// NDC → Screen
	Vector3 screenPos;
	screenPos.x = (clipPos.x + 1.0f) * 0.5f * screenWidth;
	screenPos.y = (1.0f - clipPos.y) * 0.5f * screenHeight;
	screenPos.z = clipPos.z;

	return screenPos;
}

Matrix4x4 MakeRotateAxisAngle(const Vector3& axis, float angle)
{
	Matrix4x4 matrix = {};
	float cos = std::cos(angle);
	float sin = std::sin(angle);
	float t = 1.0f - cos;
	matrix.m[0][0] = t * axis.x * axis.x + cos;
	matrix.m[0][1] = t * axis.x * axis.y + sin * axis.z;
	matrix.m[0][2] = t * axis.x * axis.z - sin * axis.y;
	matrix.m[0][3] = 0.0f;
	matrix.m[1][0] = t * axis.x * axis.y - sin * axis.z;
	matrix.m[1][1] = t * axis.y * axis.y + cos;
	matrix.m[1][2] = t * axis.y * axis.z + sin * axis.x;
	matrix.m[1][3] = 0.0f;
	matrix.m[2][0] = t * axis.x * axis.z + sin * axis.y;
	matrix.m[2][1] = t * axis.y * axis.z - sin * axis.x;
	matrix.m[2][2] = t * axis.z * axis.z + cos;
	matrix.m[2][3] = 0.0f;
	matrix.m[3][0] = 0.0f;
	matrix.m[3][1] = 0.0f;
	matrix.m[3][2] = 0.0f;
	matrix.m[3][3] = 1.0f;
	return matrix;
}

Matrix4x4 DirectionToDirection(const Vector3& from, const Vector3& to)
{
	Vector3 u = Normalize(from);
	Vector3 v = Normalize(to);

	float cosTheta = Dot(u, v);

	if (cosTheta > 0.9999f)
	{
		return MakeIdentityMatrix();
	}

	if (cosTheta < -0.9999f)
	{
		Vector3 axis = { u.y, -u.x, 0 };
	
		if (std::fabs(axis.x) < 0.0001f && std::fabs(axis.y) < 0.0001f) 
		{
			axis = Normalize(Cross(u, Vector3{ 1, 0, 0 }));

			if (std::fabs(axis.x) < 0.0001f && std::fabs(axis.y) < 0.0001f && std::fabs(axis.z) < 0.0001f) {
				axis = Normalize(Cross(u, Vector3{ 0, 1, 0 }));
			}
		}
		else
		{
			axis = Normalize(axis); 
		}

		return MakeRotateAxisAngle(axis, 3.14159265f);
	}

	Vector3 axis = Normalize(Cross(u, v));
	float angle = acosf(cosTheta);

	return MakeRotateAxisAngle(axis, angle);
}

Quaternion MakeRotateAxisAngleQuaternion(const Vector3& axis, float angle)
{
	// 軸を正規化する
	Vector3 normalizedAxis = Normalize(axis);

	// 角度を半分にする
	float halfAngle = angle * 0.5f;

	// sinとcosの値を計算
	float sinHalfAngle = sinf(halfAngle);
	float cosHalfAngle = cosf(halfAngle);

	Quaternion result;
	result.w = cosHalfAngle;
	result.x = normalizedAxis.x * sinHalfAngle;
	result.y = normalizedAxis.y * sinHalfAngle;
	result.z = normalizedAxis.z * sinHalfAngle;

	return result;
}

Vector3 RotateVector(const Vector3& vector, const Quaternion& quaternion)
{
	// ベクトルを純粋クォータニオンとして表現
	Quaternion p = { vector.x, vector.y, vector.z, 0.0f };

	Quaternion conj = Conjugate(quaternion); 

	Quaternion temp = Multiply(quaternion, p);
	Quaternion rotated_p = Multiply(temp, conj);

	// 結果のクォータニオンのベクトル部分が回転後のベクトル
	Vector3 result;
	result.x = rotated_p.x;
	result.y = rotated_p.y;
	result.z = rotated_p.z;

	return result;
}

Matrix4x4 MakeRotateMatrixQuaternion(const Quaternion& quaternion)
{
	float x = quaternion.x;
	float y = quaternion.y;
	float z = quaternion.z;
	float w = quaternion.w;

	Matrix4x4 result;

	// 1行目 
	result.m[0][0] = 1.0f - 2.0f * y * y - 2.0f * z * z;
	result.m[0][1] = 2.0f * x * y + 2.0f * w * z;
	result.m[0][2] = 2.0f * x * z - 2.0f * w * y;
	result.m[0][3] = 0.0f;                       

	// 2行目 
	result.m[1][0] = 2.0f * x * y - 2.0f * w * z;
	result.m[1][1] = 1.0f - 2.0f * x * x - 2.0f * z * z;
	result.m[1][2] = 2.0f * y * z + 2.0f * w * x;
	result.m[1][3] = 0.0f;                       

	// 3行目 
	result.m[2][0] = 2.0f * x * z + 2.0f * w * y;
	result.m[2][1] = 2.0f * y * z - 2.0f * w * x;
	result.m[2][2] = 1.0f - 2.0f * x * x - 2.0f * y * y;
	result.m[2][3] = 0.0f;                       

	// 4行目 
	result.m[3][0] = 0.0f;
	result.m[3][1] = 0.0f;
	result.m[3][2] = 0.0f;
	result.m[3][3] = 1.0f;

	return result;
}

float Dot(const Quaternion& q1, const Quaternion& q2) {
	return q1.w * q2.w + q1.x * q2.x + q1.y * q2.y + q1.z * q2.z;
}

Quaternion Slerp(const Quaternion& q0, const Quaternion& q1, float t)
{
	// tを0.0fから1.0fの範囲にクランプ
	if (t < 0.0f) t = 0.0f;
	if (t > 1.0f) t = 1.0f;

	float cosTheta = Dot(q0, q1);
	Quaternion finalQ1 = q1; 

	// 最短経路で補間するために、ドット積が負の場合にq1の符号を反転させる
	if (cosTheta < 0.0f)
	{
		cosTheta = -cosTheta;
		finalQ1.w = -finalQ1.w;
		finalQ1.x = -finalQ1.x;
		finalQ1.y = -finalQ1.y;
		finalQ1.z = -finalQ1.z;
	}

	// ほとんど同じクォータニオンの場合、Lerp（線形補間）にフォールバック
	const float THRESHOLD = 0.9995f; 
	if (cosTheta > THRESHOLD)
	{
		// 線形補間 (Lerp)
		Quaternion result;
		result.w = q0.w * (1.0f - t) + finalQ1.w * t;
		result.x = q0.x * (1.0f - t) + finalQ1.x * t;
		result.y = q0.y * (1.0f - t) + finalQ1.y * t;
		result.z = q0.z * (1.0f - t) + finalQ1.z * t;
		return Normalize(result);
	}

	// 角度を計算
	float theta = acosf(cosTheta);

	// Slerpの計算
	float sinTheta = sinf(theta);
	float invSinTheta = 1.0f / sinTheta;

	float coeff0 = sinf((1.0f - t) * theta) * invSinTheta;
	float coeff1 = sinf(t * theta) * invSinTheta;

	Quaternion result;
	result.w = q0.w * coeff0 + finalQ1.w * coeff1;
	result.x = q0.x * coeff0 + finalQ1.x * coeff1;
	result.y = q0.y * coeff0 + finalQ1.y * coeff1;
	result.z = q0.z * coeff0 + finalQ1.z * coeff1;

	return result;
}