﻿#include "Matrix.h"
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

	// q * p * q_inv を計算 (q_inv は単位クォータニオンの場合 Conjugate(q) )
	Quaternion q_inv = Conjugate(quaternion); // 単位クォータニオンであると仮定

	Quaternion temp = Multiply(quaternion, p);
	Quaternion rotated_p = Multiply(temp, q_inv);

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

	float x2 = x + x;
	float y2 = y + y;
	float z2 = z + z;
	float xx = x * x2;
	float xy = x * y2;
	float xz = x * z2;
	float yy = y * y2;
	float yz = y * z2;
	float zz = z * z2;
	float wx = w * x2;
	float wy = w * y2;
	float wz = w * z2;

	Matrix4x4 result = MakeIdentityMatrix(); // 単位行列で初期化

	// 1行目
	result.m[0][0] = 1.0f - (yy + zz);
	result.m[0][1] = xy - wz;
	result.m[0][2] = xz + wy;
	result.m[0][3] = 0.0f; // 回転行列なので、通常最後の列は [0 0 0 1]^T の形

	// 2行目
	result.m[1][0] = xy + wz;
	result.m[1][1] = 1.0f - (xx + zz);
	result.m[1][2] = yz - wx;
	result.m[1][3] = 0.0f;

	// 3行目
	result.m[2][0] = xz - wy;
	result.m[2][1] = yz + wx;
	result.m[2][2] = 1.0f - (xx + yy);
	result.m[2][3] = 0.0f;

	// 4行目
	result.m[3][0] = 0.0f;
	result.m[3][1] = 0.0f;
	result.m[3][2] = 0.0f;
	result.m[3][3] = 1.0f;

	return result;
}