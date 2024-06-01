#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <unordered_map>
#include <cmath>
#include <cctype>

using namespace std;
using namespace Eigen;

// MIT License
//
// Copyright (c) 2024 John W. Melkumov
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

//$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ DimerGen #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// Author: John W. Melkumov (jmelkumov@gmail.com).                          #
// Date: 06.01.2024.                                                        #
// Description:                                                             #
// This program takes as input two XYZ files (See Ref. * below.)            #
// with each file containing atom labels and coordinates of one of          #
// two monomers, A and B, as well as the desired center-of-mass to          #
// center-of-mass separation between them in angstroms, and 5 Euler         #
// angles in the ZYZ convention.                                            #
// Note: Although 6 Euler angles are used, (alphaA, betaA, gammaA) and      #
// (alphaB, betaB, gammaB), alphaA is hardcoded and set to 0.               #
// Ref.:                                                                    #
// * www.ccl.net/chemistry/resources/messages/1996/10/21.005-dir/index.html #
//$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$#$

unordered_map<string, float> massdict = {
    {"O", 16.000}, {"H", 1.000}, {"C", 12.000}, {"N", 14.000}, 
    {"M", 0.000}, {"L", 0.000}, {"He", 4.003}, {"Li", 6.940},
    {"Be", 9.012}, {"B", 10.810}, {"F", 18.998}, {"Ne", 20.180},
    {"Na", 22.990}, {"Mg", 24.305}, {"Al", 26.982}, {"Si", 28.085},
    {"P", 30.974}, {"S", 32.060}, {"Cl", 35.450}, {"Ar", 39.948},
    {"K", 39.098}, {"Ca", 40.078}, {"Sc", 44.956}, {"Ti", 47.867},
    {"V", 50.942}, {"Cr", 51.996}, {"Mn", 54.938}, {"Fe", 55.845},
    {"Co", 58.933}, {"Ni", 58.693}, {"Cu", 63.546}, {"Zn", 65.380},
    {"Ga", 69.723}, {"Ge", 72.630}, {"As", 74.922}, {"Se", 78.971},
    {"Br", 79.904}, {"Kr", 83.798}, {"Rb", 85.468}, {"Sr", 87.620},
    {"Y", 88.906}, {"Zr", 91.224}, {"Nb", 92.906}, {"Mo", 95.950},
    {"Tc", 98.000}, {"Ru", 101.070}, {"Rh", 102.906}, {"Pd", 106.420},
    {"Ag", 107.868}, {"Cd", 112.414}, {"In", 114.818}, {"Sn", 118.710},
    {"Sb", 121.760}, {"Te", 127.600}, {"I", 126.904}, {"Xe", 131.293},
    {"Cs", 132.905}, {"Ba", 137.327}, {"La", 138.905}, {"Ce", 140.116},
    {"Pr", 140.907}, {"Nd", 144.242}, {"Pm", 145.000}, {"Sm", 150.360},
    {"Eu", 151.964}, {"Gd", 157.250}, {"Tb", 158.925}, {"Dy", 162.500},
    {"Ho", 164.930}, {"Er", 167.259}, {"Tm", 168.934}, {"Yb", 173.045},
    {"Lu", 174.966}, {"Hf", 178.490}, {"Ta", 180.948}, {"W", 183.840},
    {"Re", 186.207}, {"Os", 190.230}, {"Ir", 192.217}, {"Pt", 195.084},
    {"Au", 196.967}, {"Hg", 200.592}, {"Tl", 204.383}, {"Pb", 207.200},
    {"Bi", 208.980}, {"Po", 209.000}, {"At", 210.000}, {"Rn", 222.000},
    {"Fr", 223.000}, {"Ra", 226.000}, {"Ac", 227.000}, {"Th", 232.038},
    {"Pa", 231.036}, {"U", 238.029}, {"Np", 237.000}, {"Pu", 244.000},
    {"Am", 243.000}, {"Cm", 247.000}, {"Bk", 247.000}, {"Cf", 251.000},
    {"Es", 252.000}, {"Fm", 257.000}, {"Md", 258.000}, {"No", 259.000},
    {"Lr", 266.000}, {"Rf", 267.000}, {"Db", 270.000}, {"Sg", 271.000},
    {"Bh", 270.000}, {"Hs", 269.000}, {"Mt", 278.000}, {"Ds", 281.000},
    {"Rg", 282.000}, {"Cn", 285.000}, {"Nh", 286.000}, {"Fl", 289.000},
    {"Mc", 290.000}, {"Lv", 293.000}, {"Ts", 294.000}, {"Og", 294.000}
};

vector<float> getMasses(const vector<string>& labels) {
    vector<float> masses;
    for (const auto& label : labels) {
        if (label.size() > 1 && islower(label[1])) {
            masses.push_back(massdict[label.substr(0, 2)]);
        } else {
            masses.push_back(massdict[label.substr(0, 1)]);
        }
    }
    return masses;
}

void readXYZ(const string& filename, vector<string>& labels, MatrixXf& xyz) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        exit(EXIT_FAILURE);
    }

    string line;
    getline(file, line); // Skip the first line of XYZ file
    getline(file, line); // Skip the second line XYZ file

    vector<vector<float>> temp_xyz;
    while (getline(file, line)) {
        if (line.length() >= 4) {
            istringstream iss(line);
            string label;
            float x, y, z;
            iss >> label >> x >> y >> z;
            labels.push_back(label);
            temp_xyz.push_back({x, y, z});
        }
    }

    xyz = MatrixXf(temp_xyz.size(), 3);
    for (size_t i = 0; i < temp_xyz.size(); ++i) {
        xyz(i, 0) = temp_xyz[i][0];
        xyz(i, 1) = temp_xyz[i][1];
        xyz(i, 2) = temp_xyz[i][2];
    }
}

Vector3f computeCOM(const MatrixXf& xyz, const vector<float>& masses) {
    Vector3f com(0.0f, 0.0f, 0.0f);
    float total_mass = 0.0f;

    for (int i = 0; i < xyz.rows(); ++i) {
        com += xyz.row(i).transpose() * masses[i];
        total_mass += masses[i];
    }

    return com / total_mass;
}

Matrix3f computeInertiaTensor(const MatrixXf& xyz, const vector<float>& masses) {
    Matrix3f inertia = Matrix3f::Zero();

    for (int i = 0; i < xyz.rows(); ++i) {
        Vector3f r = xyz.row(i).transpose();
        inertia += masses[i] * (r.dot(r) * Matrix3f::Identity() - r * r.transpose());
    }

    return inertia;
}

Matrix3f getRotationMatrix(float alpha, float beta, float gamma) {
    AngleAxisf firstZRotation(alpha * M_PI / 180.0, Vector3f::UnitZ());
    AngleAxisf yRotation(beta * M_PI / 180.0, Vector3f::UnitY());
    AngleAxisf secondZRotation(gamma * M_PI / 180.0, Vector3f::UnitZ());

    Quaternion<float> q = secondZRotation * yRotation * firstZRotation;
    return q.matrix();
}

void sortEigenvaluesAndVectors(Vector3f& eigenvalues, Matrix3f& eigenvectors) {
    vector<pair<float, Vector3f>> eigenPairs;
    for (int i = 0; i < 3; ++i) {
        eigenPairs.push_back(make_pair(eigenvalues[i], eigenvectors.col(i)));
    }

    sort(eigenPairs.begin(), eigenPairs.end(),
         [](const pair<float, Vector3f>& a, const pair<float, Vector3f>& b) {
             return a.first < b.first;
         });

    for (int i = 0; i < 3; ++i) {
        eigenvalues[i] = eigenPairs[i].first;
        eigenvectors.col(i) = eigenPairs[i].second;
    }
}

int main(int argc, char* argv[]) {
    if (argc < 9) {
        cout << "Usage: " << argv[0] << " xyzfileA xyzfileB sep betaA gammaA alphaB betaB gammaB" << endl;
        return EXIT_FAILURE;
    }

    string xyzfileA = argv[1];
    string xyzfileB = argv[2];
    float sep = stof(argv[3]);
    float alphaA = 0.0;
    float betaA = stof(argv[4]);
    float gammaA = stof(argv[5]);
    float alphaB = stof(argv[6]);
    float betaB = stof(argv[7]);
    float gammaB = stof(argv[8]);

    vector<string> labelsA, labelsB;
    MatrixXf xyzA, xyzB;

    readXYZ(xyzfileA, labelsA, xyzA);
    readXYZ(xyzfileB, labelsB, xyzB);

    vector<float> massA = getMasses(labelsA);
    vector<float> massB = getMasses(labelsB);

    Vector3f comA = computeCOM(xyzA, massA);
    Vector3f comB = computeCOM(xyzB, massB);

    MatrixXf centered_monomer_A = xyzA.rowwise() - comA.transpose();
    MatrixXf centered_monomer_B = xyzB.rowwise() - comB.transpose();

    Matrix3f inertia_A = computeInertiaTensor(centered_monomer_A, massA);
    Matrix3f inertia_B = computeInertiaTensor(centered_monomer_B, massB);

    SelfAdjointEigenSolver<Matrix3f> solverA(inertia_A);
    SelfAdjointEigenSolver<Matrix3f> solverB(inertia_B);

    Vector3f eigenvaluesA = solverA.eigenvalues();
    Matrix3f eigenvectorsA = solverA.eigenvectors();
    sortEigenvaluesAndVectors(eigenvaluesA, eigenvectorsA);

    Vector3f eigenvaluesB = solverB.eigenvalues();
    Matrix3f eigenvectorsB = solverB.eigenvectors();
    sortEigenvaluesAndVectors(eigenvaluesB, eigenvectorsB);

    Matrix3f principal_axes_A = eigenvectorsA;
    Matrix3f principal_axes_B = eigenvectorsB;

    MatrixXf monA_principal = centered_monomer_A * principal_axes_A.transpose();
    MatrixXf monB_principal = centered_monomer_B * principal_axes_B.transpose();

    Matrix3f rotation_matA = getRotationMatrix(alphaA, betaA, gammaA);
    Matrix3f rotation_matB = getRotationMatrix(alphaB, betaB, gammaB);

    MatrixXf monA_rotated_principal = monA_principal * rotation_matA;
    MatrixXf monB_rotated_principal = monB_principal * rotation_matB;

    MatrixXf centered_monA_rotated = monA_rotated_principal * principal_axes_A;
    MatrixXf centered_monB_rotated = monB_rotated_principal * principal_axes_B;

    Vector3f RAB = comB - comA;
    Vector3f RABuvec = RAB.normalized();
    Vector3f translation = sep * RABuvec;

    MatrixXf trans_rot_monomer_B = centered_monB_rotated.rowwise() + translation.transpose();

    cout << centered_monA_rotated.rows() + trans_rot_monomer_B.rows() << endl;
    cout << "Dimer_" << sep << "_Å_COM-COM_separation_" << alphaA << "deg_" << betaA << "deg_" << gammaA << "deg_" << alphaB << "deg_" << betaB << "deg_" << gammaB << "deg" << endl;

    for (int i = 0; i < centered_monA_rotated.rows(); ++i) {
        cout << labelsA[i] << " "
             << fixed << setprecision(6)
             << centered_monA_rotated(i, 0) << " "
             << centered_monA_rotated(i, 1) << " "
             << centered_monA_rotated(i, 2) << endl;
    }

    for (int i = 0; i < trans_rot_monomer_B.rows(); ++i) {
        cout << labelsB[i] << " "
             << fixed << setprecision(6)
             << trans_rot_monomer_B(i, 0) << " "
             << trans_rot_monomer_B(i, 1) << " "
             << trans_rot_monomer_B(i, 2) << endl;
    }

    return 0;
}
