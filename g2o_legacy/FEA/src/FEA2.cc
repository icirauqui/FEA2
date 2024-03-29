/**
 * This file is an addon to ORB-SLAM2.
 *
 * Its goal is to create a mesh from the 3D points within the computed map.
 * Also, it will assembly a rigidity matrix for such mesh, and compute the
 * strain deformation energy.
 *
 * Copyright (C) 2017 Íñigo Cirauqui Viloria <icirauqui at gmail dot com>
 * (University of Zaragoza)
 *
 * This addon, as ORB-SLAM2, is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * This addon is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with ORB-SLAM2. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * This file is part of ORB-SLAM2.
 *
 * Copyright (C) 2014-2016 Raúl Mur-Artal <raulmur at unizar dot es> (University
 * of Zaragoza) For more information see <https://github.com/raulmur/ORB_SLAM2>
 *
 * ORB-SLAM2 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ORB-SLAM2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with ORB-SLAM2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../include/FEA2.h"

FEA2::FEA2(unsigned int input_nFrameId, unsigned int input_E, float input_nu,
           float input_h, float input_fg1, int input_nElType, bool bSetDebug)
    : nFrameId(input_nFrameId),
      E(input_E),
      nu(input_nu),
      h(input_h),
      fg(input_fg1),
      nElType(input_nElType),
      bDebugMode(bSetDebug) {
  Ksize = 0.0;
  sE = 0.0;

  lambda = (nu * E) / ((1 + nu) * (1 - 2 * nu));
  G = E / (2 * (1 + nu));

  D.clear();
  D.push_back(vector<float>({(lambda + 2 * G), lambda, lambda, 0.0, 0.0, 0.0}));
  D.push_back(vector<float>({lambda, (lambda + 2 * G), lambda, 0.0, 0.0, 0.0}));
  D.push_back(vector<float>({lambda, lambda, (lambda + 2 * G), 0.0, 0.0, 0.0}));
  D.push_back(vector<float>({0.0, 0.0, 0.0, G, 0.0, 0.0}));
  D.push_back(vector<float>({0.0, 0.0, 0.0, 0.0, G, 0.0}));
  D.push_back(vector<float>({0.0, 0.0, 0.0, 0.0, 0.0, G}));

  gs.clear();
  gs.push_back(vector<float>({-fg, -fg, -fg}));
  gs.push_back(vector<float>({+fg, -fg, -fg}));
  gs.push_back(vector<float>({+fg, +fg, -fg}));
  gs.push_back(vector<float>({-fg, +fg, -fg}));
  gs.push_back(vector<float>({-fg, -fg, +fg}));
  gs.push_back(vector<float>({+fg, -fg, +fg}));
  gs.push_back(vector<float>({+fg, +fg, +fg}));
  gs.push_back(vector<float>({-fg, +fg, +fg}));

  vpMPs_pairs_ = std::vector<std::vector<int>>(2, std::vector<int>());
  qvecs_ = std::vector<Eigen::Vector4d>(2);
  tvecs_ = std::vector<Eigen::Vector3d>(2);
}

FEA2::~FEA2() {}

bool FEA2::Compute() {
  GetMapPointCoordinates();

  // Load points into cloud
  pcl::PointCloud<pcl::PointXYZ>::Ptr pc_0 =
      std::make_shared<pcl::PointCloud<pcl::PointXYZ>>();

  pc_0->width = vMPsXYZN_t.size();  // cloudTracked.width    = 5;
  pc_0->height = 1;
  pc_0->is_dense = false;
  pc_0->points.resize(pc_0->width * pc_0->height);

  for (size_t i = 0; i < pc_0->points.size(); i++) {
    pc_0->points[i].x = vMPsXYZN_t[i][0];
    pc_0->points[i].y = vMPsXYZN_t[i][1];
    pc_0->points[i].z = vMPsXYZN_t[i][2];
  }

  // MLS

  // pcl::PointCloud<pcl::PointXYZ>::ConstPtr pc_0p =
  // std::unique_ptr<pcl::PointCloud<PointXYZ>>(&pc_0);
  // pcl::PointCloud<pcl::PointXYZ>::ConstPtr pc_0p (new
  // pcl::PointCloud<pcl::PointXYZ> (pc_0)); return false;
  pcl::search::KdTree<pcl::PointXYZ>::Ptr tree1(
      new pcl::search::KdTree<pcl::PointXYZ>);
  pcl::MovingLeastSquares<pcl::PointXYZ, pcl::PointNormal> mls;

  // Create KD-Tree & MLS & set parameters
  mls.setComputeNormals(true);
  mls.setInputCloud(pc_0);
  // mls.setPolynomialFit (true);
  mls.setSearchMethod(tree1);
  mls.setSearchRadius(mls_search_radius_);
  mls.setPolynomialOrder(mls_polynomial_order_);

  pcl::PointCloud<pcl::PointNormal> pc_1;
  mls.process(pc_1);

  // Get corresponding indexes: for each output point, returns the index of the
  // input one.
  PointIndicesPtr pIdx1 = mls.getCorrespondingIndices();
  for (unsigned int i = 0; i < pIdx1->indices.size(); i++)
    vtindices.push_back(pIdx1->indices[i]);

  int idxit = 0;
  for (unsigned int i = 0; i < vMPsXYZN_t.size(); i++) {
    int currentpos = i;
    if (currentpos == vtindices[idxit]) {
      vbtMPsActive[i] = true;
      idxit++;
    } else
      vbtMPsActive[i] = false;
  }

  // delete pc_0p;
  // delete tree1;

  // pc_0p = nullptr;
  // tree1 = nullptr;

  // Compute Mesh triangulating by Greedy Projection
  PointCloud<PointNormal>::Ptr pc_1p(new pcl::PointCloud<PointNormal>(pc_1));

  // Search tree
  pcl::search::KdTree<pcl::PointNormal>::Ptr tree2(
      new pcl::search::KdTree<pcl::PointNormal>);
  tree2->setInputCloud(pc_1p);

  // Greedy Projection
  static pcl::GreedyProjectionTriangulation<pcl::PointNormal> GreedyProj3;

  GreedyProj3.setMu(mesh_mu_);
  GreedyProj3.setSearchRadius(mesh_search_radius_);
  GreedyProj3.setMaximumNearestNeighbors(mesh_max_neighbours_);
  GreedyProj3.setMaximumSurfaceAngle(mesh_surf_angle_ * M_PI / 180);
  GreedyProj3.setMinimumAngle(mesh_min_angle_ * M_PI / 180);
  GreedyProj3.setMaximumAngle(mesh_max_angle_ * M_PI / 180);
  GreedyProj3.setNormalConsistency(true);

  GreedyProj3.setInputCloud(pc_1p);
  GreedyProj3.setSearchMethod(tree2);

  pcl::PolygonMesh mesh;
  GreedyProj3.reconstruct(mesh);

  if (mesh.polygons.size() > 0) {
    // Get the vertex indexes of each triangle
    for (unsigned int i = 0; i < mesh.polygons.size(); i++) {
      unsigned int nver0 = mesh.polygons[i].vertices[0];
      unsigned int nver1 = mesh.polygons[i].vertices[1];
      unsigned int nver2 = mesh.polygons[i].vertices[2];

      vector<int> triangle;
      triangle.push_back(nver0);
      triangle.push_back(nver1);
      triangle.push_back(nver2);
      triangles_t.push_back(triangle);

      vector<Point3D *> vpMP2Draw;
      vpMP2Draw.push_back(vpMPs_t[vtindices[nver0]]);
      vpMP2Draw.push_back(vpMPs_t[vtindices[nver1]]);
      vpMP2Draw.push_back(vpMPs_t[vtindices[nver2]]);
      vpMPs2Draw.push_back(vpMP2Draw);
    }
  } else {
    cout << "No data" << endl;
    return false;
  }

  // if (nElType==1)         //C3D6
  //     tri2quad();

  SetSecondLayer(1);
  Set_u0(vpMPs_t, 1);
  // MatrixAssemblyC3D8(1);
  MatrixAssemblyC3D6(1);
  ImposeDirichletEncastre_K(1, vvDir_t, 100000000.0);

  return true;
}

void FEA2::ComputeTestMesh() {
  GetMapPointCoordinates();

  if (LoadMPsIntoCloud())
    if (MLS()) ComputeMesh();
}

bool FEA2::Compute(int nMode) {
  GetMapPointCoordinates();

  if (LoadMPsIntoCloud(nMode)) {
    if (bDebugMode) cout << "           - MLS" << endl;
    if (MLS(nMode)) {
      if (bDebugMode) cout << "           - Tri mesh" << endl;
      if (ComputeMesh(nMode)) {
        if (nElType == 1)  // C3D6
          tri2quad();

        SetSecondLayer(nMode);

        if (nMode == 1) Set_u0(vpMPs_t, nMode);

        if (nElType == 1) MatrixAssemblyC3D8(nMode);
        if (nElType == 2) MatrixAssemblyC3D6(nMode);

        if (nMode == 1) ImposeDirichletEncastre_K(nMode, vvDir_t, 100000000.0);

        return true;
      }
    }
  }
  return false;
}

bool FEA2::ComputePlanar() {
  if (LoadMPsIntoCloud()) {
    if (Triangulate()) {
      if (nElType == 1)  // C3D8
        tri2quad();

      SetSecondLayer();

      Set_u0(vpMPs_t);

      if (nElType == 1) MatrixAssemblyC3D8();
      if (nElType == 2) MatrixAssemblyC3D6();

      ImposeDirichletEncastre_K(1, vvDir_t, 100000000.0);

      return true;
    }
  }

  return false;
}

void FEA2::GetMapPointCoordinates() {
  for (unsigned int i = 0; i < vpMPs_t.size(); i++) {
    Point3D *pMP = vpMPs_t[i];
    Eigen::Vector3d vPos = pMP->XYZ();

    vector<float> vfMP;
    vfMP.push_back(vPos(0));
    vfMP.push_back(vPos(1));
    vfMP.push_back(vPos(2));

    vMPsXYZN_t.push_back(vfMP);
    vbtMPsActive.push_back(true);
  }
}

void FEA2::ClearMapPointCoordinates() {
  vMPsXYZN_t.clear();
  vbtMPsActive.clear();
  vpMPs2Draw.clear();
  // triangles_t.clear();
}

void FEA2::AddMapPointCoordinates(std::vector<float> vPoint) {
  vMPsXYZN_t.push_back(vPoint);
  vbtMPsActive.push_back(true);
}

bool FEA2::LoadMPsIntoCloud(int nMode) {
  // Load top layer into cloud
  pc_t_0.width = vMPsXYZN_t.size();  // cloudTracked.width    = 5;
  pc_t_0.height = 1;
  pc_t_0.is_dense = false;
  pc_t_0.points.resize(pc_t_0.width * pc_t_0.height);

  for (size_t i = 0; i < pc_t_0.points.size(); i++) {
    pc_t_0.points[i].x = vMPsXYZN_t[i][0];
    pc_t_0.points[i].y = vMPsXYZN_t[i][1];
    pc_t_0.points[i].z = vMPsXYZN_t[i][2];
  }

  if (pc_t_0.width > 0)
    return true;
  else
    return false;
}

bool FEA2::MLS(int nMode) {
  // pc_0 = pcl::PointCloud<pcl::PointXYZ>::Ptr(new pcl::PointCloud<PointXYZ>
  // (pc_t_0));
  pcl::PointCloud<pcl::PointXYZ>::Ptr pc_01(
      new pcl::PointCloud<PointXYZ>(pc_t_0));

  // Create KD-Tree & MLS & set parameters
  pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(
      new pcl::search::KdTree<pcl::PointXYZ>);
  pcl::MovingLeastSquares<pcl::PointXYZ, pcl::PointNormal> mls;
  mls.setComputeNormals(true);
  mls.setInputCloud(pc_01);
  // mls.setPolynomialFit (true);
  mls.setSearchMethod(tree);
  mls.setSearchRadius(mls_search_radius_);  // 900000 //100 //0.25    //
                                            // Original 0.03 // Default 0.00
  mls.setPolynomialOrder(mls_polynomial_order_);  // 3 //50 // Default 2

  mls.process(pc_t_1);

  // Get corresponding indexes: for each output point, returns the index of the
  // input one.
  PointIndicesPtr pIdx1 = mls.getCorrespondingIndices();
  for (unsigned int i = 0; i < pIdx1->indices.size(); i++)
    vtindices.push_back(pIdx1->indices[i]);

  int idxit = 0;
  for (unsigned int i = 0; i < vMPsXYZN_t.size(); i++) {
    int currentpos = i;
    if (currentpos == vtindices[idxit]) {
      vbtMPsActive[i] = true;
      idxit++;
    } else
      vbtMPsActive[i] = false;
  }

  if (pc_t_1.width > 0)
    return true;
  else
    return false;
}

bool FEA2::ComputeMesh(int nMode) {
  // TRIANGULATION, GREEDY PROJECTION
  PointCloud<PointNormal>::Ptr ptr_pc_t_mesh_1(
      new pcl::PointCloud<PointNormal>(pc_t_1));

  // Search tree
  pcl::search::KdTree<pcl::PointNormal>::Ptr tree(
      new pcl::search::KdTree<pcl::PointNormal>);
  tree->setInputCloud(ptr_pc_t_mesh_1);

  // Greedy Projection parameters
  int nInMu = 5;
  int nSearchRad = 30;
  int nMaxNeig = 70;
  int nSurAng = 150;
  int nMinAng = 45;
  int nMaxAng = 90;

  // CalculateGP3Parameters(ptr_pc_t_mesh_1,&nInMu,&nSearchRad,&nMaxNeig,&nSurAng,&nMinAng,&nMaxAng);
  // if (bDebugMode) cout << "           -
  // GP3Parameters(InMu/nSearchRad/nMaxNeig/nSurAng/nMinAng/nMaxAng) = "
  //                      << nInMu << " / " << nSearchRad << " / " << nMaxNeig
  //                      << " / " << nSurAng << " / " << nMinAng << " / " <<
  //                      nMaxAng << endl;

  // Greedy Projection
  static pcl::GreedyProjectionTriangulation<pcl::PointNormal> GreedyProj3;

  GreedyProj3.setMu(mesh_mu_);                       // original 2.8
  GreedyProj3.setSearchRadius(mesh_search_radius_);  // original 0.50
  GreedyProj3.setMaximumNearestNeighbors(mesh_max_neighbours_);  // original 150
  GreedyProj3.setMaximumSurfaceAngle(mesh_surf_angle_ * M_PI /
                                     180);                    // 45 degrees
  GreedyProj3.setMinimumAngle(mesh_min_angle_ * M_PI / 180);  // 10 degrees
  GreedyProj3.setMaximumAngle(mesh_max_angle_ * M_PI / 180);  // 120 degrees 2/3
  GreedyProj3.setNormalConsistency(true);

  GreedyProj3.setInputCloud(ptr_pc_t_mesh_1);
  GreedyProj3.setSearchMethod(tree);
  GreedyProj3.reconstruct(mesh_t);

  if (mesh_t.polygons.size() > 0) {
    // Get the vertex indexes of each triangle
    for (unsigned int i = 0; i < mesh_t.polygons.size(); i++) {
      unsigned int nver0 = mesh_t.polygons[i].vertices[0];
      unsigned int nver1 = mesh_t.polygons[i].vertices[1];
      unsigned int nver2 = mesh_t.polygons[i].vertices[2];

      vector<int> triangle;
      triangle.push_back(nver0);
      triangle.push_back(nver1);
      triangle.push_back(nver2);
      triangles_t.push_back(triangle);

      vector<Point3D *> vpMP2Draw;
      vpMP2Draw.push_back(vpMPs_t[vtindices[nver0]]);
      vpMP2Draw.push_back(vpMPs_t[vtindices[nver1]]);
      vpMP2Draw.push_back(vpMPs_t[vtindices[nver2]]);
      vpMPs2Draw.push_back(vpMP2Draw);
    }

    return true;
  } else {
    cout << "No data" << endl;
    return false;
  }
}

bool FEA2::Triangulate() {
  // Cloud must be filled before calling the function
  if (pc_t_0.width == 0) {
    return false;
  }

  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(
      new pcl::PointCloud<pcl::PointXYZ>(pc_t_0));

  // Normal estimation*
  pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> n;
  pcl::PointCloud<pcl::Normal>::Ptr normals(new pcl::PointCloud<pcl::Normal>);
  pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(
      new pcl::search::KdTree<pcl::PointXYZ>);
  tree->setInputCloud(cloud);
  n.setInputCloud(cloud);
  n.setSearchMethod(tree);
  n.setKSearch(20);
  n.compute(*normals);
  //* normals should not contain the point normals + surface curvatures

  // Concatenate the XYZ and normal fields*
  pcl::PointCloud<pcl::PointNormal>::Ptr cloud_with_normals(
      new pcl::PointCloud<pcl::PointNormal>);
  pcl::concatenateFields(*cloud, *normals, *cloud_with_normals);
  //* cloud_with_normals = cloud + normals

  // Create search tree*
  pcl::search::KdTree<pcl::PointNormal>::Ptr tree2(
      new pcl::search::KdTree<pcl::PointNormal>);
  tree2->setInputCloud(cloud_with_normals);

  // Initialize objects
  pcl::GreedyProjectionTriangulation<pcl::PointNormal> gp3;
  pcl::PolygonMesh mesh_out;

  // Set the maximum distance between connected points (maximum edge length)
  gp3.setSearchRadius(500.0);

  // Set typical values for the parameters
  gp3.setMu(2.5);
  gp3.setMaximumNearestNeighbors(100);
  gp3.setMaximumSurfaceAngle(M_PI / 4);  // 45 degrees
  gp3.setMinimumAngle(M_PI / 18);        // 10 degrees
  gp3.setMaximumAngle(2 * M_PI / 3);     // 120 degrees
  gp3.setNormalConsistency(false);

  // Get result
  gp3.setInputCloud(cloud_with_normals);
  gp3.setSearchMethod(tree2);
  gp3.reconstruct(mesh_out);

  // Additional vertex information
  // To which part (cloud) does each vertex belong
  std::vector<int> parts = gp3.getPartIDs();
  // Whether the vertex status is [-1,0,1,2,3] =
  // [NONE,FREE,FRINGE,BOUNDARY,COMPLETED]
  std::vector<int> states = gp3.getPointStates();

  // Get largest part
  std::unordered_map<int, int> map;
  int max_val = -1;
  int max_key = -1;
  for (unsigned int i = 0; i < parts.size(); i++) {
    if (map.find(parts[i]) == map.end()) {
      map[parts[i]] = 1;
      if (max_val < 1) {
        max_val = 1;
        max_key = parts[i];
      }
    } else {
      map[parts[i]]++;
      if (map[parts[i]] > max_val) {
        max_val = map[parts[i]];
        max_key = parts[i];
      }
    }
  }

  // Get triangles of largest part only
  triangles_t.clear();
  std::vector<std::vector<int>> triangles;
  for (unsigned int i = 0; i < mesh_out.polygons.size(); i++) {
    unsigned int nver0 = mesh_out.polygons[i].vertices[0];
    unsigned int nver1 = mesh_out.polygons[i].vertices[1];
    unsigned int nver2 = mesh_out.polygons[i].vertices[2];

    if (parts[nver0] == max_key && parts[nver1] == max_key &&
        parts[nver2] == max_key) {
      std::vector<int> triangle;
      triangle.push_back(nver0);
      triangle.push_back(nver1);
      triangle.push_back(nver2);
      triangles_t.push_back(triangle);
    }
  }

  if (triangles_t.size() > 0)
    return true;
  else
    return false;
}

std::vector<std::vector<int>> FEA2::GetMeshTriangles(int nMode) {
  return triangles_t;
}

void FEA2::CalculateGP3Parameters(pcl::PointCloud<pcl::PointNormal>::Ptr ppc,
                                  int *pnInMu, int *pnSearchRad, int *pnMaxNeig,
                                  int *pnSurAng, int *pnMinAng, int *pnMaxAng) {
  vector<vector<float>> vit1;
  for (unsigned int i = 0; i < ppc->points.size(); i++) {
    vector<float> viti;
    viti.push_back(ppc->points[i].x);
    viti.push_back(ppc->points[i].y);
    viti.push_back(ppc->points[i].z);
    vit1.push_back(viti);
  }

  vector<float> vtot, vdmin, vdmax;
  float dd, dx, dy, dz;
  for (unsigned int i = 0; i < vit1.size(); i++) {
    float dmin = 1000.0;
    float dmax = 0.0;
    for (unsigned int j = 0; j < vit1.size(); j++) {
      if (j == i) continue;

      dx = vit1[j][0] - vit1[i][0];
      dy = vit1[j][1] - vit1[i][1];
      dz = vit1[j][2] - vit1[i][2];
      dx *= dx;
      dy *= dy;
      dz *= dz;
      dd = dx + dy + dz;
      dd = sqrt(dd);

      if (dd < dmin) dmin = dd;
      if (dd > dmax) dmax = dd;

      vtot.push_back(dd);
    }
    vdmin.push_back(dmin);
    vdmax.push_back(dmax);
  }

  float fmedtot = 0.0;
  float fsigmatot = 0.0;
  for (unsigned int i = 0; i < vtot.size(); i++) {
    fmedtot += vtot[i] / vtot.size();
  }
  for (unsigned int i = 0; i < vtot.size(); i++) {
    fsigmatot += (vtot[i] - fmedtot) * (vtot[i] - fmedtot);
  }
  fsigmatot /= vtot.size();
  fsigmatot = sqrt(fsigmatot);

  float ftminmin = 1000.0;
  float ftmaxmin = 0.0;
  float ftminmax = 1000.0;
  float ftmaxmax = 0.0;
  for (unsigned int i = 0; i < vdmin.size(); i++) {
    if (vdmin[i] < ftminmin) ftminmin = vdmin[i];
    if (vdmin[i] > ftmaxmin) ftmaxmin = vdmin[i];
    if (vdmax[i] < ftminmax) ftminmax = vdmax[i];
    if (vdmax[i] > ftmaxmax) ftmaxmax = vdmax[i];
  }

  float fmedmin = (ftmaxmin - ftminmin) / 2;
  float fmedmax = (ftmaxmax - ftminmax) / 2;
  float fsigmamin = 0.0;
  float fsigmamax = 0.0;
  for (unsigned int i = 0; i < vdmin.size(); i++) {
    fsigmamin += (vdmin[i] - fmedmin) * (vdmin[i] - fmedmin);
    fsigmamax += (vdmax[i] - fmedmax) * (vdmax[i] - fmedmax);
    fsigmamin /= vdmin.size();
    fsigmamax /= vdmax.size();
    fsigmamin = sqrt(fsigmamin);
    fsigmamax = sqrt(fsigmamax);
  }

  *pnInMu = ceil(fmedtot / fmedmin);
  *pnSearchRad = ceil(fmedtot);
  *pnMaxNeig = vtot.size() / 4;  // fmedtot + fsigmatot;
  *pnSurAng = 150;
  *pnMinAng = 45;
  *pnMaxAng = 180 - 2 * 45;
}

void FEA2::tri2quad() {
  // Get all vertices in mesh
  vector<int> vertices;
  for (unsigned int i = 0; i < triangles_t.size(); i++)
    for (unsigned int j = 0; j < triangles_t[i].size(); j++)
      if (!(std::find(vertices.begin(), vertices.end(), triangles_t[i][j]) !=
            vertices.end()))
        vertices.push_back(triangles_t[i][j]);

  int maxvertex = vertices.back();

  vector<int> offvertices;
  for (int i = 0; i < maxvertex; i++)
    if (!(std::find(vertices.begin(), vertices.end(), i) != vertices.end()))
      offvertices.push_back(i);

  for (unsigned int i = 0; i < vertices.size(); i++)
    vbtMPsActive[vtindices[vertices[i]]] = true;

  for (unsigned int i = 0; i < offvertices.size(); i++)
    vbtMPsActive[vtindices[offvertices[i]]] = false;

  // Get true indexes
  vector<vector<int>> tempidx;
  for (unsigned int i = 0; i < triangles_t.size(); i++) {
    vector<int> triangle;
    for (unsigned int j = 0; j < triangles_t[i].size(); j++)
      triangle.push_back(vtindices[triangles_t[i][j]]);
    tempidx.push_back(triangle);
  }

  triangles_t.clear();
  for (unsigned int i = 0; i < tempidx.size(); i++)
    triangles_t.push_back(tempidx[i]);

  // Record the 3 edges_t in each triangle
  for (unsigned int i = 0; i < triangles_t.size(); i++) {
    if (triangles_t[i].size() < 3) continue;

    vector<int> edge1;
    vector<int> edge2;
    vector<int> edge3;
    edge1.push_back(triangles_t[i][0]);
    edge1.push_back(triangles_t[i][1]);
    edge2.push_back(triangles_t[i][0]);
    edge2.push_back(triangles_t[i][2]);
    edge3.push_back(triangles_t[i][1]);
    edge3.push_back(triangles_t[i][2]);

    vector<vector<int>> tempedges;
    tempedges.push_back(edge1);
    tempedges.push_back(edge2);
    tempedges.push_back(edge3);

    edgespertriangle_t.push_back(tempedges);
  }

  vector<vector<vector<int>>> vCommonEdges;
  vCommonEdges.resize(triangles_t.size());
  for (unsigned int i = 0; i < vCommonEdges.size(); i++)
    vCommonEdges[i].resize(3);
  /*triangulo1
          commonedge1
                  edge1 triangulo2 edge2
          commonedge2
                  edge1 triangulo2 edge2
          commonedge3
                  edge1 triangulo2 edge2
  triangulo2
  triangulo3
  ...*/

  vIdxNewVertices.clear();
  vIdxNewVertices.resize(triangles_t.size());
  // for(unsigned int i=0; i<vIdxNewVertices.size(); i++)
  // for (unsigned int j=0; j<4; j++)
  // vIdxNewVertices[i].push_back(-1);
  /*triangulo1
          idxvertexedge0 idxvertexedge1 idxvertexedge2 idxbaricentro
  triangulo2
  triangulo3
  ...*/

  vector<bool> vbComputedTriangle = vector<bool>(triangles_t.size(), false);
  /*triangulo1
  triangulo2
  triangulo3
  ...*/

  quads_t.clear();

  // Fill vCommonEdges with the function originally designed to join common
  // edges.
  for (unsigned int i = 0; i < edgespertriangle_t.size(); i++) {
    // Won't happen, but just in case
    if (edgespertriangle_t[i].size() < 3) continue;

    // For each triangle, we look for common edges with all the other polygons
    for (unsigned int j = 0; j < edgespertriangle_t.size(); j++) {
      // Don't look within itself
      if (j == i) continue;

      // Again, won't happen, but just in case
      if (edgespertriangle_t[j].size() < 3) continue;

      // For all the edges within the reference triangle
      for (unsigned int k = 0; k < edgespertriangle_t[i].size(); k++) {
        // Compare the first edge in the reference triangle with all the
        // three(k) edges in the candidate
        if (((edgespertriangle_t[i][0][0] == edgespertriangle_t[j][k][0]) &&
             (edgespertriangle_t[i][0][1] == edgespertriangle_t[j][k][1])) ||
            ((edgespertriangle_t[i][0][0] == edgespertriangle_t[j][k][1]) &&
             (edgespertriangle_t[i][0][1] == edgespertriangle_t[j][k][0]))) {
          // edge1 = 0 & 1
          vCommonEdges[i][0].push_back(0);
          vCommonEdges[i][0].push_back(j);
          vCommonEdges[i][0].push_back(k);
        }

        // Compare the second edge in the reference triangle with all the
        // three(k) edges in the candidate
        if (((edgespertriangle_t[i][1][0] == edgespertriangle_t[j][k][0]) &&
             (edgespertriangle_t[i][1][1] == edgespertriangle_t[j][k][1])) ||
            ((edgespertriangle_t[i][1][0] == edgespertriangle_t[j][k][1]) &&
             (edgespertriangle_t[i][1][1] == edgespertriangle_t[j][k][0]))) {
          // edge1 = 0 & 2
          vCommonEdges[i][1].push_back(1);
          vCommonEdges[i][1].push_back(j);
          vCommonEdges[i][1].push_back(k);
        }

        // Compare the third edge in the reference triangle with all the
        // three(k) edges in the candidate
        if (((edgespertriangle_t[i][2][0] == edgespertriangle_t[j][k][0]) &&
             (edgespertriangle_t[i][2][1] == edgespertriangle_t[j][k][1])) ||
            ((edgespertriangle_t[i][2][0] == edgespertriangle_t[j][k][1]) &&
             (edgespertriangle_t[i][2][1] == edgespertriangle_t[j][k][0]))) {
          // edge1 = 1 & 2
          vCommonEdges[i][2].push_back(2);
          vCommonEdges[i][2].push_back(j);
          vCommonEdges[i][2].push_back(k);
        }
        // k changes here, so in each loop, we compare all the edges in
        // the reference triangle with only one edge of the candidate
      }
    }  // End of the loop through all the candidates for a reference triangle
  }    // End of the loop through all the triangles set as reference

  vNewPointsBase.clear();

  for (unsigned int i = 0; i < triangles_t.size(); i++) {
    // Get vertices and indices for refferece triangle
    int v0idx = triangles_t[i][0];
    vector<float> v0 = vMPsXYZN_t[v0idx];

    int v1idx = triangles_t[i][1];
    vector<float> v1 = vMPsXYZN_t[v1idx];

    int v2idx = triangles_t[i][2];
    vector<float> v2 = vMPsXYZN_t[v2idx];

    // Set variables to store the ressult
    vector<float> m01;
    int m01idx;
    vector<float> m02;
    int m02idx;
    vector<float> m12;
    int m12idx;
    vector<float> g;
    int gidx;
    vector<int> quad1;
    vector<int> quad2;
    vector<int> quad3;

    // Compute edge 0
    for (unsigned int j = 0; j < vCommonEdges[i].size(); j++) {
      vector<int> NewPointsBasei;
      if (!vCommonEdges[i][j].empty()) {
        if (vCommonEdges[i][j][0] == 0) {
          if (vbComputedTriangle[vCommonEdges[i][j][1]]) {
            m01idx =
                vIdxNewVertices[vCommonEdges[i][j][1]][vCommonEdges[i][j][2]];
            m01 = vMPsXYZN_t[m01idx];
          } else {
            m01idx = vMPsXYZN_t.size();

            m01.push_back((v0[0] + v1[0]) / 2);
            m01.push_back((v0[1] + v1[1]) / 2);
            m01.push_back((v0[2] + v1[2]) / 2);

            vMPsXYZN_t.push_back(m01);

            NewPointsBasei.push_back(v0idx);
            NewPointsBasei.push_back(v1idx);
            vNewPointsBase.push_back(NewPointsBasei);
          }
        }
      } else {
        m01idx = vMPsXYZN_t.size();

        m01.push_back((v0[0] + v1[0]) / 2);
        m01.push_back((v0[1] + v1[1]) / 2);
        m01.push_back((v0[2] + v1[2]) / 2);

        vMPsXYZN_t.push_back(m01);

        NewPointsBasei.push_back(v0idx);
        NewPointsBasei.push_back(v1idx);
        vNewPointsBase.push_back(NewPointsBasei);
      }
    }

    // Compute edge 1
    for (unsigned int j = 0; j < vCommonEdges[i].size(); j++) {
      vector<int> NewPointsBasei;
      if (!vCommonEdges[i][j].empty()) {
        if (vCommonEdges[i][j][0] == 1) {
          if (vbComputedTriangle[vCommonEdges[i][j][1]]) {
            m02idx =
                vIdxNewVertices[vCommonEdges[i][j][1]][vCommonEdges[i][j][2]];
            m02 = vMPsXYZN_t[m02idx];
          } else {
            m02idx = vMPsXYZN_t.size();

            m02.push_back((v0[0] + v2[0]) / 2);
            m02.push_back((v0[1] + v2[1]) / 2);
            m02.push_back((v0[2] + v2[2]) / 2);

            vMPsXYZN_t.push_back(m02);

            NewPointsBasei.push_back(v0idx);
            NewPointsBasei.push_back(v2idx);
            vNewPointsBase.push_back(NewPointsBasei);
          }
        }
      } else {
        m02idx = vMPsXYZN_t.size();

        m02.push_back((v0[0] + v2[0]) / 2);
        m02.push_back((v0[1] + v2[1]) / 2);
        m02.push_back((v0[2] + v2[2]) / 2);

        vMPsXYZN_t.push_back(m02);

        NewPointsBasei.push_back(v0idx);
        NewPointsBasei.push_back(v2idx);
        vNewPointsBase.push_back(NewPointsBasei);
      }
    }

    // Compute edge 2
    for (unsigned int j = 0; j < vCommonEdges[i].size(); j++) {
      vector<int> NewPointsBasei;
      if (!vCommonEdges[i][j].empty()) {
        if (vCommonEdges[i][j][0] == 2) {
          if (vbComputedTriangle[vCommonEdges[i][j][1]]) {
            m12idx =
                vIdxNewVertices[vCommonEdges[i][j][1]][vCommonEdges[i][j][2]];
            m12 = vMPsXYZN_t[m12idx];
          } else {
            m12idx = vMPsXYZN_t.size();

            m12.push_back((v1[0] + v2[0]) / 2);
            m12.push_back((v1[1] + v2[1]) / 2);
            m12.push_back((v1[2] + v2[2]) / 2);

            vMPsXYZN_t.push_back(m12);

            NewPointsBasei.push_back(v1idx);
            NewPointsBasei.push_back(v2idx);
            vNewPointsBase.push_back(NewPointsBasei);
          }
        }
      } else {
        m12idx = vMPsXYZN_t.size();

        m12.push_back((v1[0] + v2[0]) / 2);
        m12.push_back((v1[1] + v2[1]) / 2);
        m12.push_back((v1[2] + v2[2]) / 2);

        vMPsXYZN_t.push_back(m12);

        NewPointsBasei.push_back(v1idx);
        NewPointsBasei.push_back(v2idx);
        vNewPointsBase.push_back(NewPointsBasei);
      }
    }

    // Compute barycentre in each triangle
    g.push_back((v0[0] + v1[0] + v2[0]) / 3);
    g.push_back((v0[1] + v1[1] + v2[1]) / 3);
    g.push_back((v0[2] + v1[2] + v2[2]) / 3);
    gidx = vMPsXYZN_t.size();
    vMPsXYZN_t.push_back(g);

    vector<int> NewPointsBasei;
    NewPointsBasei.push_back(v0idx);
    NewPointsBasei.push_back(v1idx);
    NewPointsBasei.push_back(v2idx);
    vNewPointsBase.push_back(NewPointsBasei);

    // Set the three quads
    quad1.push_back(v0idx);
    quad1.push_back(m01idx);
    quad1.push_back(gidx);
    quad1.push_back(m02idx);

    quad2.push_back(v1idx);
    quad2.push_back(m12idx);
    quad2.push_back(gidx);
    quad2.push_back(m01idx);

    quad3.push_back(v2idx);
    quad3.push_back(m02idx);
    quad3.push_back(gidx);
    quad3.push_back(m12idx);

    // Add the quads to the global vector
    quads_t.push_back(quad1);
    quads_t.push_back(quad2);
    quads_t.push_back(quad3);

    // Set vertices idx
    vIdxNewVertices[i].push_back(m01idx);
    vIdxNewVertices[i].push_back(m02idx);
    vIdxNewVertices[i].push_back(m12idx);
    vIdxNewVertices[i].push_back(gidx);

    // Set triangle as already computed
    vbComputedTriangle[i] = true;
  }
}

void FEA2::SetSecondLayer(int nMode) {
  vvDir_t.clear();
  for (unsigned int i = 0; i < vMPsXYZN_t.size(); i++) {
    vector<float> point;
    float pos1 = vMPsXYZN_t[i][0] - h;
    float pos2 = vMPsXYZN_t[i][1] - h;
    float pos3 = vMPsXYZN_t[i][2] - h;
    point.push_back(pos1);
    point.push_back(pos2);
    point.push_back(pos3);
    vMPsXYZN_t2.push_back(point);

    vector<int> vDiri;
    vDiri.push_back(vMPsXYZN_t.size() + i);
    vvDir_t.push_back(vDiri);
  }
}

void FEA2::Set_u0(vector<Point3D *> vpMPs, int nMode) {
  u0.clear();
  for (unsigned int i = 0; i < vMPsXYZN_t.size(); i++)
    for (unsigned int j = 0; j < vMPsXYZN_t[i].size(); j++)
      u0.push_back(vMPsXYZN_t[i][j]);
  for (unsigned int i = 0; i < vMPsXYZN_t2.size(); i++)
    for (unsigned int j = 0; j < vMPsXYZN_t2[i].size(); j++)
      u0.push_back(vMPsXYZN_t2[i][j]);
}

vector<vector<float>> FEA2::ComputeKeiC3D8(vector<vector<float>> vfPts) {
  vector<vector<float>> vBtDB =
      vector<vector<float>>(24, vector<float>(24, 0.0));

  for (unsigned int ops = 0; ops < gs.size(); ops++) {
    float xi = gs[ops][0];
    float eta = gs[ops][1];
    float zeta = gs[ops][2];

    float dN1dxi = -0.125 * ((1 - eta) * (1 - zeta));
    float dN1deta = -0.125 * ((1 - xi) * (1 - zeta));
    float dN1dzeta = -0.125 * ((1 - xi) * (1 - eta));
    float dN2dxi = +0.125 * ((1 - eta) * (1 - zeta));
    float dN2deta = -0.125 * ((1 + xi) * (1 - zeta));
    float dN2dzeta = -0.125 * ((1 + xi) * (1 - eta));
    float dN3dxi = +0.125 * ((1 + eta) * (1 - zeta));
    float dN3deta = +0.125 * ((1 + xi) * (1 - zeta));
    float dN3dzeta = -0.125 * ((1 + xi) * (1 + eta));
    float dN4dxi = -0.125 * ((1 + eta) * (1 - zeta));
    float dN4deta = +0.125 * ((1 - xi) * (1 - zeta));
    float dN4dzeta = -0.125 * ((1 - xi) * (1 + eta));
    float dN5dxi = -0.125 * ((1 - eta) * (1 + zeta));
    float dN5deta = -0.125 * ((1 - xi) * (1 + zeta));
    float dN5dzeta = +0.125 * ((1 - xi) * (1 - eta));
    float dN6dxi = +0.125 * ((1 - eta) * (1 + zeta));
    float dN6deta = -0.125 * ((1 + xi) * (1 + zeta));
    float dN6dzeta = +0.125 * ((1 + xi) * (1 - eta));
    float dN7dxi = +0.125 * ((1 + eta) * (1 + zeta));
    float dN7deta = +0.125 * ((1 + xi) * (1 + zeta));
    float dN7dzeta = +0.125 * ((1 + xi) * (1 + eta));
    float dN8dxi = -0.125 * ((1 + eta) * (1 + zeta));
    float dN8deta = +0.125 * ((1 - xi) * (1 + zeta));
    float dN8dzeta = +0.125 * ((1 - xi) * (1 + eta));

    /*dxdxi*/ float J_00 = dN1dxi * vfPts[0][0] + dN2dxi * vfPts[1][0] +
                           dN3dxi * vfPts[2][0] + dN4dxi * vfPts[3][0] +
                           dN5dxi * vfPts[4][0] + dN6dxi * vfPts[5][0] +
                           dN7dxi * vfPts[6][0] + dN8dxi * vfPts[7][0];
    /*dydxi*/ float J_01 = dN1dxi * vfPts[0][1] + dN2dxi * vfPts[1][1] +
                           dN3dxi * vfPts[2][1] + dN4dxi * vfPts[3][1] +
                           dN5dxi * vfPts[4][1] + dN6dxi * vfPts[5][1] +
                           dN7dxi * vfPts[6][1] + dN8dxi * vfPts[7][1];
    /*dzdxi*/ float J_02 = dN1dxi * vfPts[0][2] + dN2dxi * vfPts[1][2] +
                           dN3dxi * vfPts[2][2] + dN4dxi * vfPts[3][2] +
                           dN5dxi * vfPts[4][2] + dN6dxi * vfPts[5][2] +
                           dN7dxi * vfPts[6][2] + dN8dxi * vfPts[7][2];
    /*dxdeta*/ float J_10 = dN1deta * vfPts[0][0] + dN2deta * vfPts[1][0] +
                            dN3deta * vfPts[2][0] + dN4deta * vfPts[3][0] +
                            dN5deta * vfPts[4][0] + dN6deta * vfPts[5][0] +
                            dN7deta * vfPts[6][0] + dN8deta * vfPts[7][0];
    /*dydeta*/ float J_11 = dN1deta * vfPts[0][1] + dN2deta * vfPts[1][1] +
                            dN3deta * vfPts[2][1] + dN4deta * vfPts[3][1] +
                            dN5deta * vfPts[4][1] + dN6deta * vfPts[5][1] +
                            dN7deta * vfPts[6][1] + dN8deta * vfPts[7][1];
    /*dzdeta*/ float J_12 = dN1deta * vfPts[0][2] + dN2deta * vfPts[1][2] +
                            dN3deta * vfPts[2][2] + dN4deta * vfPts[3][2] +
                            dN5deta * vfPts[4][2] + dN6deta * vfPts[5][2] +
                            dN7deta * vfPts[6][2] + dN8deta * vfPts[7][2];
    /*dxdzeta*/ float J_20 = dN1dzeta * vfPts[0][0] + dN2dzeta * vfPts[1][0] +
                             dN3dzeta * vfPts[2][0] + dN4dzeta * vfPts[3][0] +
                             dN5dzeta * vfPts[4][0] + dN6dzeta * vfPts[5][0] +
                             dN7dzeta * vfPts[6][0] + dN8dzeta * vfPts[7][0];
    /*dydzeta*/ float J_21 = dN1dzeta * vfPts[0][1] + dN2dzeta * vfPts[1][1] +
                             dN3dzeta * vfPts[2][1] + dN4dzeta * vfPts[3][1] +
                             dN5dzeta * vfPts[4][1] + dN6dzeta * vfPts[5][1] +
                             dN7dzeta * vfPts[6][1] + dN8dzeta * vfPts[7][1];
    /*dzdzeta*/ float J_22 = dN1dzeta * vfPts[0][2] + dN2dzeta * vfPts[1][2] +
                             dN3dzeta * vfPts[2][2] + dN4dzeta * vfPts[3][2] +
                             dN5dzeta * vfPts[4][2] + dN6dzeta * vfPts[5][2] +
                             dN7dzeta * vfPts[6][2] + dN8dzeta * vfPts[7][2];

    float Jac = J_00 * J_11 * J_22 + J_01 * J_12 * J_20 + J_10 * J_21 * J_02 -
                J_20 * J_11 * J_02 - J_10 * J_01 * J_22 - J_21 * J_12 * J_00;

    float J1_00 = (+1) * ((J_11 * J_22) - (J_21 * J_12)) / Jac;
    float J1_01 = (-1) * ((J_01 * J_22) - (J_21 * J_02)) / Jac;
    float J1_02 = (-1) * ((J_01 * J_12) - (J_11 * J_02)) / Jac;
    float J1_10 = (-1) * ((J_10 * J_22) - (J_20 * J_12)) / Jac;
    float J1_11 = (-1) * ((J_00 * J_22) - (J_20 * J_02)) / Jac;
    float J1_12 = (-1) * ((J_00 * J_12) - (J_10 * J_02)) / Jac;
    float J1_20 = (+1) * ((J_10 * J_21) - (J_20 * J_11)) / Jac;
    float J1_21 = (-1) * ((J_00 * J_21) - (J_20 * J_01)) / Jac;
    float J1_22 = (-1) * ((J_00 * J_11) - (J_10 * J_01)) / Jac;

    float dN1dx = J1_00 * dN1dxi + J1_01 * dN1deta + J1_02 * dN1dzeta;
    float dN1dy = J1_10 * dN1dxi + J1_11 * dN1deta + J1_12 * dN1dzeta;
    float dN1dz = J1_20 * dN1dxi + J1_21 * dN1deta + J1_22 * dN1dzeta;
    float dN2dx = J1_00 * dN2dxi + J1_01 * dN2deta + J1_02 * dN2dzeta;
    float dN2dy = J1_10 * dN2dxi + J1_11 * dN2deta + J1_12 * dN2dzeta;
    float dN2dz = J1_20 * dN2dxi + J1_21 * dN2deta + J1_22 * dN2dzeta;
    float dN3dx = J1_00 * dN3dxi + J1_01 * dN3deta + J1_02 * dN3dzeta;
    float dN3dy = J1_10 * dN3dxi + J1_11 * dN3deta + J1_12 * dN3dzeta;
    float dN3dz = J1_20 * dN3dxi + J1_21 * dN3deta + J1_22 * dN3dzeta;
    float dN4dx = J1_00 * dN4dxi + J1_01 * dN4deta + J1_02 * dN4dzeta;
    float dN4dy = J1_10 * dN4dxi + J1_11 * dN4deta + J1_12 * dN4dzeta;
    float dN4dz = J1_20 * dN4dxi + J1_21 * dN4deta + J1_22 * dN4dzeta;
    float dN5dx = J1_00 * dN5dxi + J1_01 * dN5deta + J1_02 * dN5dzeta;
    float dN5dy = J1_10 * dN5dxi + J1_11 * dN5deta + J1_12 * dN5dzeta;
    float dN5dz = J1_20 * dN5dxi + J1_21 * dN5deta + J1_22 * dN5dzeta;
    float dN6dx = J1_00 * dN6dxi + J1_01 * dN6deta + J1_02 * dN6dzeta;
    float dN6dy = J1_10 * dN6dxi + J1_11 * dN6deta + J1_12 * dN6dzeta;
    float dN6dz = J1_20 * dN6dxi + J1_21 * dN6deta + J1_22 * dN6dzeta;
    float dN7dx = J1_00 * dN7dxi + J1_01 * dN7deta + J1_02 * dN7dzeta;
    float dN7dy = J1_10 * dN7dxi + J1_11 * dN7deta + J1_12 * dN7dzeta;
    float dN7dz = J1_20 * dN7dxi + J1_21 * dN7deta + J1_22 * dN7dzeta;
    float dN8dx = J1_00 * dN8dxi + J1_01 * dN8deta + J1_02 * dN8dzeta;
    float dN8dy = J1_10 * dN8dxi + J1_11 * dN8deta + J1_12 * dN8dzeta;
    float dN8dz = J1_20 * dN8dxi + J1_21 * dN8deta + J1_22 * dN8dzeta;

    float B[6][24] = {
        dN1dx, 0.0,   0.0,   dN2dx, 0.0,   0.0,   dN3dx, 0.0,   0.0,   dN4dx,
        0.0,   0.0,   dN5dx, 0.0,   0.0,   dN6dx, 0.0,   0.0,   dN7dx, 0.0,
        0.0,   dN8dx, 0.0,   0.0,   0.0,   dN1dy, 0.0,   0.0,   dN2dy, 0.0,
        0.0,   dN3dy, 0.0,   0.0,   dN4dy, 0.0,   0.0,   dN5dy, 0.0,   0.0,
        dN6dy, 0.0,   0.0,   dN7dy, 0.0,   0.0,   dN8dy, 0.0,   0.0,   0.0,
        dN1dz, 0.0,   0.0,   dN2dz, 0.0,   0.0,   dN3dz, 0.0,   0.0,   dN4dz,
        0.0,   0.0,   dN5dz, 0.0,   0.0,   dN6dz, 0.0,   0.0,   dN7dz, 0.0,
        0.0,   dN8dz, dN1dy, dN1dx, 0.0,   dN2dy, dN2dx, 0.0,   dN3dy, dN3dx,
        0.0,   dN4dy, dN4dx, 0.0,   dN5dy, dN5dx, 0.0,   dN6dy, dN6dx, 0.0,
        dN7dy, dN7dx, 0.0,   dN8dy, dN8dx, 0.0,   dN1dz, 0.0,   dN1dx, dN2dz,
        0.0,   dN2dx, dN3dz, 0.0,   dN3dx, dN4dz, 0.0,   dN4dx, dN5dz, 0.0,
        dN5dx, dN6dz, 0.0,   dN6dx, dN7dz, 0.0,   dN7dx, dN8dz, 0.0,   dN8dx,
        0.0,   dN1dz, dN1dy, 0.0,   dN2dz, dN2dy, 0.0,   dN3dz, dN3dy, 0.0,
        dN4dz, dN4dy, 0.0,   dN5dz, dN5dy, 0.0,   dN6dz, dN6dy, 0.0,   dN7dz,
        dN7dy, 0.0,   dN8dz, dN8dy};

    float BtD[24][6] = {0.0};

    for (unsigned int i = 0; i < 24; i++)
      for (unsigned int j = 0; j < 6; j++)
        BtD[i][j] = B[0][i] * D[0][j] + B[1][i] * D[1][j] + B[2][i] * D[2][j] +
                    B[3][i] * D[3][j] + B[4][i] * D[4][j] + B[5][i] * D[5][j];

    for (unsigned int i = 0; i < 24; i++)
      for (unsigned int j = 0; j < 24; j++) {
        float vBtDBaux = BtD[i][0] * B[0][j] + BtD[i][1] * B[1][j] +
                         BtD[i][2] * B[2][j] + BtD[i][3] * B[3][j] +
                         BtD[i][4] * B[4][j] + BtD[i][5] * B[5][j];
        vBtDB[i][j] += vBtDBaux * Jac;
      }
  }

  return vBtDB;
}

vector<vector<float>> FEA2::ComputeKeiC3D6(vector<vector<float>> vfPts) {
  vector<vector<float>> vBtDB =
      vector<vector<float>>(18, vector<float>(18, 0.0));

  for (unsigned int ops = 0; ops < gs.size(); ops++) {
    float xi = gs[ops][0];
    float eta = gs[ops][1];
    float zeta = gs[ops][2];

    float dN1dxi = -(1 + zeta) / 2;
    float dN1deta = -(1 + zeta) / 2;
    float dN1dzeta = (1 - xi - eta) / 2;
    float dN2dxi = (1 + zeta) / 2;
    float dN2deta = 0.0;
    float dN2dzeta = xi / 2;
    float dN3dxi = 0.0;
    float dN3deta = (1 + zeta) / 2;
    float dN3dzeta = eta / 2;
    float dN4dxi = -(1 - zeta) / 2;
    float dN4deta = -(1 - zeta) / 2;
    float dN4dzeta = -(1 - xi - eta) / 2;
    float dN5dxi = (1 - zeta) / 2;
    float dN5deta = 0.0;
    float dN5dzeta = -xi / 2;
    float dN6dxi = 0.0;
    float dN6deta = (1 - zeta) / 2;
    float dN6dzeta = -eta / 2;

    /*dxdxi*/ float J_00 = dN1dxi * vfPts[0][0] + dN2dxi * vfPts[1][0] +
                           dN3dxi * vfPts[2][0] + dN4dxi * vfPts[3][0] +
                           dN5dxi * vfPts[4][0] + dN6dxi * vfPts[5][0];
    /*dydxi*/ float J_01 = dN1dxi * vfPts[0][1] + dN2dxi * vfPts[1][1] +
                           dN3dxi * vfPts[2][1] + dN4dxi * vfPts[3][1] +
                           dN5dxi * vfPts[4][1] + dN6dxi * vfPts[5][1];
    /*dzdxi*/ float J_02 = dN1dxi * vfPts[0][2] + dN2dxi * vfPts[1][2] +
                           dN3dxi * vfPts[2][2] + dN4dxi * vfPts[3][2] +
                           dN5dxi * vfPts[4][2] + dN6dxi * vfPts[5][2];
    /*dxdeta*/ float J_10 = dN1deta * vfPts[0][0] + dN2deta * vfPts[1][0] +
                            dN3deta * vfPts[2][0] + dN4deta * vfPts[3][0] +
                            dN5deta * vfPts[4][0] + dN6deta * vfPts[5][0];
    /*dydeta*/ float J_11 = dN1deta * vfPts[0][1] + dN2deta * vfPts[1][1] +
                            dN3deta * vfPts[2][1] + dN4deta * vfPts[3][1] +
                            dN5deta * vfPts[4][1] + dN6deta * vfPts[5][1];
    /*dzdeta*/ float J_12 = dN1deta * vfPts[0][2] + dN2deta * vfPts[1][2] +
                            dN3deta * vfPts[2][2] + dN4deta * vfPts[3][2] +
                            dN5deta * vfPts[4][2] + dN6deta * vfPts[5][2];
    /*dxdzeta*/ float J_20 = dN1dzeta * vfPts[0][0] + dN2dzeta * vfPts[1][0] +
                             dN3dzeta * vfPts[2][0] + dN4dzeta * vfPts[3][0] +
                             dN5dzeta * vfPts[4][0] + dN6dzeta * vfPts[5][0];
    /*dydzeta*/ float J_21 = dN1dzeta * vfPts[0][1] + dN2dzeta * vfPts[1][1] +
                             dN3dzeta * vfPts[2][1] + dN4dzeta * vfPts[3][1] +
                             dN5dzeta * vfPts[4][1] + dN6dzeta * vfPts[5][1];
    /*dzdzeta*/ float J_22 = dN1dzeta * vfPts[0][2] + dN2dzeta * vfPts[1][2] +
                             dN3dzeta * vfPts[2][2] + dN4dzeta * vfPts[3][2] +
                             dN5dzeta * vfPts[4][2] + dN6dzeta * vfPts[5][2];

    float Jac = J_00 * J_11 * J_22 + J_01 * J_12 * J_20 + J_10 * J_21 * J_02 -
                J_20 * J_11 * J_02 - J_10 * J_01 * J_22 - J_21 * J_12 * J_00;

    // cout << "Jac = " << Jac << endl;

    float J1_00 = (+1) * ((J_11 * J_22) - (J_21 * J_12)) / Jac;
    float J1_01 = (-1) * ((J_01 * J_22) - (J_21 * J_02)) / Jac;
    float J1_02 = (-1) * ((J_01 * J_12) - (J_11 * J_02)) / Jac;
    float J1_10 = (-1) * ((J_10 * J_22) - (J_20 * J_12)) / Jac;
    float J1_11 = (-1) * ((J_00 * J_22) - (J_20 * J_02)) / Jac;
    float J1_12 = (-1) * ((J_00 * J_12) - (J_10 * J_02)) / Jac;
    float J1_20 = (+1) * ((J_10 * J_21) - (J_20 * J_11)) / Jac;
    float J1_21 = (-1) * ((J_00 * J_21) - (J_20 * J_01)) / Jac;
    float J1_22 = (-1) * ((J_00 * J_11) - (J_10 * J_01)) / Jac;

    float dN1dx = J1_00 * dN1dxi + J1_01 * dN1deta + J1_02 * dN1dzeta;
    float dN1dy = J1_10 * dN1dxi + J1_11 * dN1deta + J1_12 * dN1dzeta;
    float dN1dz = J1_20 * dN1dxi + J1_21 * dN1deta + J1_22 * dN1dzeta;
    float dN2dx = J1_00 * dN2dxi + J1_01 * dN2deta + J1_02 * dN2dzeta;
    float dN2dy = J1_10 * dN2dxi + J1_11 * dN2deta + J1_12 * dN2dzeta;
    float dN2dz = J1_20 * dN2dxi + J1_21 * dN2deta + J1_22 * dN2dzeta;
    float dN3dx = J1_00 * dN3dxi + J1_01 * dN3deta + J1_02 * dN3dzeta;
    float dN3dy = J1_10 * dN3dxi + J1_11 * dN3deta + J1_12 * dN3dzeta;
    float dN3dz = J1_20 * dN3dxi + J1_21 * dN3deta + J1_22 * dN3dzeta;
    float dN4dx = J1_00 * dN4dxi + J1_01 * dN4deta + J1_02 * dN4dzeta;
    float dN4dy = J1_10 * dN4dxi + J1_11 * dN4deta + J1_12 * dN4dzeta;
    float dN4dz = J1_20 * dN4dxi + J1_21 * dN4deta + J1_22 * dN4dzeta;
    float dN5dx = J1_00 * dN5dxi + J1_01 * dN5deta + J1_02 * dN5dzeta;
    float dN5dy = J1_10 * dN5dxi + J1_11 * dN5deta + J1_12 * dN5dzeta;
    float dN5dz = J1_20 * dN5dxi + J1_21 * dN5deta + J1_22 * dN5dzeta;
    float dN6dx = J1_00 * dN6dxi + J1_01 * dN6deta + J1_02 * dN6dzeta;
    float dN6dy = J1_10 * dN6dxi + J1_11 * dN6deta + J1_12 * dN6dzeta;
    float dN6dz = J1_20 * dN6dxi + J1_21 * dN6deta + J1_22 * dN6dzeta;

    float B[6][18] = {
        dN1dx, 0.0,   0.0,   dN2dx, 0.0,   0.0,   dN3dx, 0.0,   0.0,   dN4dx,
        0.0,   0.0,   dN5dx, 0.0,   0.0,   dN6dx, 0.0,   0.0,   0.0,   dN1dy,
        0.0,   0.0,   dN2dy, 0.0,   0.0,   dN3dy, 0.0,   0.0,   dN4dy, 0.0,
        0.0,   dN5dy, 0.0,   0.0,   dN6dy, 0.0,   0.0,   0.0,   dN1dz, 0.0,
        0.0,   dN2dz, 0.0,   0.0,   dN3dz, 0.0,   0.0,   dN4dz, 0.0,   0.0,
        dN5dz, 0.0,   0.0,   dN6dz, dN1dy, dN1dx, 0.0,   dN2dy, dN2dx, 0.0,
        dN3dy, dN3dx, 0.0,   dN4dy, dN4dx, 0.0,   dN5dy, dN5dx, 0.0,   dN6dy,
        dN6dx, 0.0,   dN1dz, 0.0,   dN1dx, dN2dz, 0.0,   dN2dx, dN3dz, 0.0,
        dN3dx, dN4dz, 0.0,   dN4dx, dN5dz, 0.0,   dN5dx, dN6dz, 0.0,   dN6dx,
        0.0,   dN1dz, dN1dy, 0.0,   dN2dz, dN2dy, 0.0,   dN3dz, dN3dy, 0.0,
        dN4dz, dN4dy, 0.0,   dN5dz, dN5dy, 0.0,   dN6dz, dN6dy};

    float BtD[18][6] = {0.0};

    for (unsigned int i = 0; i < 18; i++)
      for (unsigned int j = 0; j < 6; j++)
        BtD[i][j] = B[0][i] * D[0][j] + B[1][i] * D[1][j] + B[2][i] * D[2][j] +
                    B[3][i] * D[3][j] + B[4][i] * D[4][j] + B[5][i] * D[5][j];

    for (unsigned int i = 0; i < 18; i++)
      for (unsigned int j = 0; j < 18; j++) {
        float vBtDBaux = BtD[i][0] * B[0][j] + BtD[i][1] * B[1][j] +
                         BtD[i][2] * B[2][j] + BtD[i][3] * B[3][j] +
                         BtD[i][4] * B[4][j] + BtD[i][5] * B[5][j];
        vBtDB[i][j] += vBtDBaux * Jac;
      }
  }

  return vBtDB;
}

bool FEA2::MatrixAssemblyC3D8(int nMode) {
  if (nMode == 1) {
    int nTotalNodes = vMPsXYZN_t.size() + vMPsXYZN_t2.size();
    Ksize = 3 * nTotalNodes;
    if (bDebugMode)
      cout << "                - MatAssembly (hex, 24 DoF per el. Curr: "
           << quads_t.size() << " quads, " << nTotalNodes
           << " nodes, Ksize = " << Ksize << ")" << endl;

    if (Ksize <= 3) return false;

    K = vector<vector<float>>(Ksize, vector<float>(Ksize, 0.0));

    for (unsigned int i = 0; i < quads_t.size(); i++) {
      vector<int> nodes;
      nodes.push_back(quads_t[i][0]);
      nodes.push_back(quads_t[i][1]);
      nodes.push_back(quads_t[i][2]);
      nodes.push_back(quads_t[i][3]);
      nodes.push_back(quads_t[i][0] + vMPsXYZN_t.size());
      nodes.push_back(quads_t[i][1] + vMPsXYZN_t.size());
      nodes.push_back(quads_t[i][2] + vMPsXYZN_t.size());
      nodes.push_back(quads_t[i][3] + vMPsXYZN_t.size());

      vector<vector<float>> vfPts;
      for (unsigned int j = 0; j < 4; j++) {
        vector<float> vfPtsi;
        vfPtsi.push_back(vMPsXYZN_t[nodes[j]][0]);
        vfPtsi.push_back(vMPsXYZN_t[nodes[j]][1]);
        vfPtsi.push_back(vMPsXYZN_t[nodes[j]][2]);
        vfPts.push_back(vfPtsi);
      }

      for (unsigned int j = 0; j < 4; j++) {
        vector<float> vfPtsi;
        vfPtsi.push_back(vMPsXYZN_t2[nodes[j]][0]);
        vfPtsi.push_back(vMPsXYZN_t2[nodes[j]][1]);
        vfPtsi.push_back(vMPsXYZN_t2[nodes[j]][2]);
        vfPts.push_back(vfPtsi);
      }

      vector<vector<float>> Kei = ComputeKeiC3D8(vfPts);

      vector<int> mn;
      for (unsigned j = 0; j < nodes.size(); j++) {
        mn.push_back(nodes[j] * 3);
      }

      for (unsigned int ni = 0; ni < nodes.size(); ni++) {
        for (unsigned int nj = 0; nj < nodes.size(); nj++) {
          for (unsigned int m = 0; m < 3; m++) {
            for (unsigned int n = 0; n < 3; n++) {
              if ((mn[ni] + m) >= Ksize || (mn[nj] + n) >= Ksize) continue;
              K[mn[ni] + m][mn[nj] + n] += Kei[3 * ni + m][3 * nj + n];
            }
          }
        }
      }
    }
    return true;
  } else if (nMode == 2) {
    int nTotalNodes = vMPsXYZN_ut.size() + vMPsXYZN_ut2.size();
    Kusize = 3 * nTotalNodes;
    if (bDebugMode)
      cout << "                - MatAssembly (hex, 24 DoF per el. Curr: "
           << quads_u.size() << " quads, " << nTotalNodes
           << " nodes, Kusize = " << Kusize << ")" << endl;

    if (Kusize <= 3) return false;

    Ku.clear();
    vector<float> ki = vector<float>(Kusize, 0.0);
    for (unsigned int i = 0; i < Kusize; i++) Ku.push_back(ki);

    for (unsigned int i = 0; i < quads_u.size(); i++) {
      vector<int> nodes;
      nodes.push_back(quads_u[i][0]);
      nodes.push_back(quads_u[i][1]);
      nodes.push_back(quads_u[i][2]);
      nodes.push_back(quads_u[i][3]);
      nodes.push_back(quads_u[i][0] + vMPsXYZN_ut.size());
      nodes.push_back(quads_u[i][1] + vMPsXYZN_ut.size());
      nodes.push_back(quads_u[i][2] + vMPsXYZN_ut.size());
      nodes.push_back(quads_u[i][3] + vMPsXYZN_ut.size());

      vector<vector<float>> vfPts;
      for (unsigned int j = 0; j < 4; j++) {
        vector<float> vfPtsi;
        vfPtsi.push_back(vMPsXYZN_ut[nodes[j]][0]);
        vfPtsi.push_back(vMPsXYZN_ut[nodes[j]][1]);
        vfPtsi.push_back(vMPsXYZN_ut[nodes[j]][2]);
        vfPts.push_back(vfPtsi);
      }

      for (unsigned int j = 0; j < 4; j++) {
        vector<float> vfPtsi;
        vfPtsi.push_back(vMPsXYZN_ut2[nodes[j]][0]);
        vfPtsi.push_back(vMPsXYZN_ut2[nodes[j]][1]);
        vfPtsi.push_back(vMPsXYZN_ut2[nodes[j]][2]);
        vfPts.push_back(vfPtsi);
      }

      vector<vector<float>> Kei = ComputeKeiC3D8(vfPts);

      vector<int> mn;
      for (unsigned j = 0; j < nodes.size(); j++) {
        mn.push_back(nodes[j] * 3);
      }

      for (unsigned int ni = 0; ni < nodes.size(); ni++) {
        for (unsigned int nj = 0; nj < nodes.size(); nj++) {
          for (unsigned int m = 0; m < 3; m++) {
            for (unsigned int n = 0; n < 3; n++) {
              if ((mn[ni] + m) >= Ksize || (mn[nj] + n) >= Ksize) continue;
              Ku[mn[ni] + m][mn[nj] + n] += Kei[3 * ni + m][3 * nj + n];
            }
          }
        }
      }
    }
    return true;
  } else {
    return false;
  }
}

bool FEA2::MatrixAssemblyC3D6(int nMode) {
  if (nMode == 1) {
    int nTotalNodes = vMPsXYZN_t.size() + vMPsXYZN_t2.size();
    Ksize = 3 * nTotalNodes;
    if (bDebugMode)
      cout << "                - MatAssembly (tri, 18 DoF per el. Curr: "
           << triangles_t.size() << " triangles, " << nTotalNodes
           << " nodes, Ksize = " << Ksize << ")" << endl;

    if (Ksize <= 3) return false;

    K = vector<vector<float>>(Ksize, vector<float>(Ksize, 0.0));

    for (unsigned int i = 0; i < triangles_t.size(); i++) {
      vector<int> nodes;
      nodes.push_back(triangles_t[i][0]);
      nodes.push_back(triangles_t[i][1]);
      nodes.push_back(triangles_t[i][2]);
      nodes.push_back(triangles_t[i][0] + vMPsXYZN_t.size());
      nodes.push_back(triangles_t[i][1] + vMPsXYZN_t.size());
      nodes.push_back(triangles_t[i][2] + vMPsXYZN_t.size());

      vector<vector<float>> vfPts;
      for (unsigned int j = 0; j < 3; j++) {
        vector<float> vfPtsi;
        vfPtsi.push_back(vMPsXYZN_t[nodes[j]][0]);
        vfPtsi.push_back(vMPsXYZN_t[nodes[j]][1]);
        vfPtsi.push_back(vMPsXYZN_t[nodes[j]][2]);
        vfPts.push_back(vfPtsi);
      }

      for (unsigned int j = 0; j < 3; j++) {
        vector<float> vfPtsi;
        vfPtsi.push_back(vMPsXYZN_t2[nodes[j]][0]);
        vfPtsi.push_back(vMPsXYZN_t2[nodes[j]][1]);
        vfPtsi.push_back(vMPsXYZN_t2[nodes[j]][2]);
        vfPts.push_back(vfPtsi);
      }

      vector<vector<float>> Kei = ComputeKeiC3D6(vfPts);

      vector<int> mn;
      for (unsigned j = 0; j < nodes.size(); j++) {
        mn.push_back(nodes[j] * 3);
      }

      for (unsigned int ni = 0; ni < nodes.size(); ni++) {
        for (unsigned int nj = 0; nj < nodes.size(); nj++) {
          for (unsigned int m = 0; m < 3; m++) {
            for (unsigned int n = 0; n < 3; n++) {
              if ((mn[ni] + m) >= Ksize || (mn[nj] + n) >= Ksize) continue;
              K[mn[ni] + m][mn[nj] + n] += Kei[3 * ni + m][3 * nj + n];
            }
          }
        }
      }
    }
    return true;
  } else if (nMode == 2) {
    int nTotalNodes = vMPsXYZN_ut.size() + vMPsXYZN_ut2.size();
    Kusize = 3 * nTotalNodes;
    if (bDebugMode)
      cout << "                - MatAssembly (tri, 18 DoF per el. Curr: "
           << triangles_u.size() << " triangles, " << nTotalNodes
           << " nodes, Kusize = " << Kusize << ")" << endl;

    if (Kusize <= 3) return false;

    Ku.clear();
    vector<float> ki = vector<float>(Kusize, 0.0);
    for (unsigned int i = 0; i < Kusize; i++) Ku.push_back(ki);

    for (unsigned int i = 0; i < triangles_u.size(); i++) {
      vector<int> nodes;
      nodes.push_back(triangles_u[i][0]);
      nodes.push_back(triangles_u[i][1]);
      nodes.push_back(triangles_u[i][2]);
      nodes.push_back(triangles_u[i][0] + vMPsXYZN_ut.size());
      nodes.push_back(triangles_u[i][1] + vMPsXYZN_ut.size());
      nodes.push_back(triangles_u[i][2] + vMPsXYZN_ut.size());

      vector<vector<float>> vfPts;
      for (unsigned int j = 0; j < 3; j++) {
        vector<float> vfPtsi;
        vfPtsi.push_back(vMPsXYZN_ut[nodes[j]][0]);
        vfPtsi.push_back(vMPsXYZN_ut[nodes[j]][1]);
        vfPtsi.push_back(vMPsXYZN_ut[nodes[j]][2]);
        vfPts.push_back(vfPtsi);
      }

      for (unsigned int j = 0; j < 3; j++) {
        vector<float> vfPtsi;
        vfPtsi.push_back(vMPsXYZN_ut2[nodes[j]][0]);
        vfPtsi.push_back(vMPsXYZN_ut2[nodes[j]][1]);
        vfPtsi.push_back(vMPsXYZN_ut2[nodes[j]][2]);
        vfPts.push_back(vfPtsi);
      }

      vector<vector<float>> Kei = ComputeKeiC3D6(vfPts);

      vector<int> mn;
      for (unsigned j = 0; j < nodes.size(); j++) {
        mn.push_back(nodes[j] * 3);
      }

      for (unsigned int ni = 0; ni < nodes.size(); ni++) {
        for (unsigned int nj = 0; nj < nodes.size(); nj++) {
          for (unsigned int m = 0; m < 3; m++) {
            for (unsigned int n = 0; n < 3; n++) {
              if ((mn[ni] + m) >= Ksize || (mn[nj] + n) >= Ksize) continue;
              Ku[mn[ni] + m][mn[nj] + n] += Kei[3 * ni + m][3 * nj + n];
            }
          }
        }
      }
    }
    return true;
  } else {
    return false;
  }
}

void FEA2::ImposeDirichletEncastre_K(int nMode, vector<vector<int>> vD,
                                     float Klarge) {
  for (unsigned int i = 0; i < vD.size(); i++) {
    int mp0 = 3 * (vD[i][0] - 1);
    int mp1 = mp0 + 1;
    int mp2 = mp0 + 2;

    if (nMode == 1) {
      K[mp0][mp0] = Klarge;
      K[mp1][mp1] = Klarge;
      K[mp2][mp2] = Klarge;
    } else if (nMode == 2) {
      Ku[mp0][mp0] = Klarge;
      Ku[mp1][mp1] = Klarge;
      Ku[mp2][mp2] = Klarge;
    }
  }
}

void FEA2::ImposeDirichletEncastre_a(vector<vector<int>> vD, float Klarge) {
  for (unsigned int i = 0; i < vD.size(); i++) {
    int mp0 = 3 * (vD[i][0] - 1);
    int mp1 = mp0 + 1;
    int mp2 = mp0 + 2;

    vva[mp0][0] = 1 / Klarge;
    vva[mp1][0] = 1 / Klarge;
    vva[mp2][0] = 1 / Klarge;
  }
}

vector<vector<float>> FEA2::InvertMatrixEigen(vector<vector<float>> m1) {
  Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> A;
  A.resize(m1.size(), m1.size());

  for (unsigned int i = 0; i < m1.size(); i++) {
    for (unsigned int j = 0; j < m1[i].size(); j++) {
      A(i, j) = m1[i][j];
    }
  }

  float detA = A.determinant();
  Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> A1 = A.inverse();

  vector<vector<float>> Ainv;
  if (detA != 0.0) {
    for (unsigned int i = 0; i < m1.size(); i++) {
      vector<float> Ainvi;
      for (unsigned int j = 0; j < m1[i].size(); j++) {
        Ainvi.push_back(A1(i, j));
      }
      Ainv.push_back(Ainvi);
    }
  } else {
    vector<float> Ainvi;
    Ainvi.push_back(0.0);
    Ainv.push_back(Ainvi);
  }

  return Ainv;
}

vector<vector<float>> FEA2::MultiplyMatricesEigen(vector<vector<float>> m1,
                                                  vector<vector<float>> m2) {
  vector<vector<float>> pr =
      vector<vector<float>>(m1.size(), vector<float>(m2[0].size(), 0.0));

  if (m1[0].size() != m2.size()) {
    pr = vector<vector<float>>(1, vector<float>(1, 0.0));
    return pr;
  }

  Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> A;
  A.resize(m1.size(), m1[0].size());
  for (unsigned int i = 0; i < m1.size(); i++) {
    for (unsigned int j = 0; j < m1[i].size(); j++) {
      A(i, j) = m1[i][j];
    }
  }

  Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> B;
  B.resize(m2.size(), m2[0].size());
  for (unsigned int i = 0; i < m2.size(); i++) {
    for (unsigned int j = 0; j < m2[i].size(); j++) {
      B(i, j) = m2[i][j];
    }
  }

  Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> C;
  C.resize(pr.size(), pr[0].size());
  C = A * B;

  for (unsigned int i = 0; i < pr.size(); i++) {
    for (unsigned int j = 0; j < pr[i].size(); j++) {
      pr[i][j] = C(i, j);
    }
  }

  return pr;
}

void FEA2::Set_uf_self() {
  std::vector<std::vector<float>> xyzs;
  for (unsigned int i=0; i<vpMPs_ut.size(); i++) {
    std::vector<float> xyz;
    xyz.push_back(vpMPs_ut[i]->X());
    xyz.push_back(vpMPs_ut[i]->Y());
    xyz.push_back(vpMPs_ut[i]->Z());
    xyzs.push_back(xyz);
  }
  Set_uf(xyzs);
}

void FEA2::Set_uf(vector<vector<float>> vPoints) {
  vMPsXYZN_t.clear();
  vMPsXYZN_t.resize(vPoints.size());

  for (unsigned int i = 0; i < vPoints.size(); i++) {
    vector<float> vfMP;

    vfMP.push_back(vPoints[i][0]);
    vfMP.push_back(vPoints[i][1]);
    vfMP.push_back(vPoints[i][2]);

    vMPsXYZN_t[i] = vfMP;
  }

  for (unsigned int i = 0; i < vNewPointsBase.size(); i++) {
    vector<float> mi;

    if (vNewPointsBase[i].size() == 2) {
      int v0idx = vNewPointsBase[i][0];
      int v1idx = vNewPointsBase[i][1];

      vector<float> v0 = vMPsXYZN_t[v0idx];
      vector<float> v1 = vMPsXYZN_t[v1idx];

      mi.push_back((v0[0] + v1[0]) / 2);
      mi.push_back((v0[1] + v1[1]) / 2);
      mi.push_back((v0[2] + v1[2]) / 2);
    }
    if (vNewPointsBase[i].size() == 3) {
      int v0idx = vNewPointsBase[i][0];
      int v1idx = vNewPointsBase[i][1];
      int v2idx = vNewPointsBase[i][2];

      vector<float> v0 = vMPsXYZN_t[v0idx];
      vector<float> v1 = vMPsXYZN_t[v1idx];
      vector<float> v2 = vMPsXYZN_t[v2idx];

      mi.push_back((v0[0] + v1[0] + v2[0]) / 3);
      mi.push_back((v0[1] + v1[1] + v2[1]) / 3);
      mi.push_back((v0[2] + v1[2] + v2[2]) / 3);
    }

    vMPsXYZN_t.push_back(mi);
  }

  uf.clear();

  for (unsigned int i = 0; i < vMPsXYZN_t.size(); i++) {
    if (vMPsXYZN_t[i].empty())
      for (unsigned int j = 0; j < 3; j++) uf.push_back(0.0);
    else
      for (unsigned int j = 0; j < 3; j++) uf.push_back(vMPsXYZN_t[i][j]);
  }

  for (unsigned int i = 0; i < vMPsXYZN_t2.size(); i++) {
    if (vMPsXYZN_t2[i].empty())
      for (unsigned int j = 0; j < 3; j++) uf.push_back(0.0);
    else
      for (unsigned int j = 0; j < 3; j++) uf.push_back(vMPsXYZN_t2[i][j]);
  }
}

void FEA2::SetDisplacement(std::vector<std::vector<float>> vDisplacement) {
  vva.clear();
  for (unsigned int i = 0; i < vDisplacement.size(); i++) {
    vector<float> va;
    va.push_back(vDisplacement[i][0]);
    va.push_back(vDisplacement[i][1]);
    va.push_back(vDisplacement[i][2]);
    vva.push_back(va);
  }
}

void FEA2::ComputeDisplacement() {
  vva.clear();
  for (unsigned int i = 0; i < u0.size(); i++) {
    vector<float> va;
    va.push_back(uf[i] - u0[i]);
    vva.push_back(va);
  }

  ImposeDirichletEncastre_a(vvDir_t, 100000000.0);
}

void FEA2::ComputeForces() {
  vvf.clear();

  /*
  double sumK = 0.0;
  for (unsigned int i=0; i<K.size(); i++){
      for (unsigned int j=0; j<K[i].size(); j++){
          sumK += K[i][j];
      }
  }

  double sumvva = 0.0;
  for (unsigned int i=0; i<vva.size(); i++){
      for (unsigned int j=0; j<vva[i].size(); j++){
          sumvva += vva[i][j];
      }
  }

  std::cout << "sumK = " << sumK << std::endl;
  std::cout << "sumvva = " << sumvva << std::endl;

  */

  vvf = MultiplyMatricesEigen(K, vva);

  vpMPs2DrawWgt = ComputeScaledDef(vvf, 100);

  vector<vector<float>> va1 = vector_resize_cols(vva, 3);
  vector<vector<float>> vf1 = vector_resize_cols(vvf, 3);

  ve_.clear();

  for (unsigned int i = 0; i < va1.size(); i++) {
    float aa = 0.0;
    float ff = 0.0;
    for (unsigned int j = 0; j < va1[i].size(); j++) {
      aa += va1[i][j] * va1[i][j];
      ff += vf1[i][j] * vf1[i][j];
    }
    aa = sqrt(aa);
    ff = sqrt(ff);
    float en = ff * aa;

    ve_.push_back(en);
  }

  /*
  vector<vector<float> > va1 = vector_resize_cols(vva,3);

  float vamin_ = 0.0;
  for (unsigned int i=0; i<va1[ifmin_].size(); i++){
      vamin_ += va1[i][j]*va1[i][j];
  }
  vamin_ = sqrt(va2);

  float vamax_ = 0.0;
  for (unsigned int i=0; i<va1[ifmax_].size(); i++){
      vamax_ += va1[i][j]*va1[i][j];
  }
  vamax_ = sqrt(va2);
  */
}

vector<float> FEA2::ComputeScaledDef(vector<vector<float>> vf,
                                     float scalerange) {
  vector<vector<float>> vf1 = vector_resize_cols(vf, 3);

  vector<float> vf2;
  for (unsigned int i = 0; i < vf1.size(); i++) {
    float vfi2 = 0.0;
    for (unsigned int j = 0; j < vf1[i].size(); j++)
      vfi2 += vf1[i][j] * vf1[i][j];
    vf2.push_back(sqrt(vfi2));
  }

  float vfmin = 0.0;
  float vfmax = 0.0;
  for (unsigned int i = 0; i < vf2.size(); i++) {
    if (i == 0) {
      vfmin = vf2[i];
      vfmax = vf2[i];
    } else if (vf2[i] < vfmin) {
      vfmin = vf2[i];
      ifmin_ = i;
    } else if (vf2[i] > vfmax) {
      vfmax = vf2[i];
      ifmax_ = i;
    }
  }

  vfmin_ = vfmin;
  vfmax_ = vfmax;

  vfmax -= vfmin;
  vector<float> vsf;
  for (unsigned int i = 0; i < vf2.size(); i++)
    vsf.push_back((vf2[i] - vfmin) / vfmax);

  vfmin = 0.0;
  vfmax = 0.0;
  vector<float> nvsf;
  for (unsigned int i = 0; i < triangles_t.size(); i++) {
    float tricolor = 0.0;
    for (unsigned int j = 0; j < triangles_t[i].size(); j++)
      tricolor += vsf[triangles_t[i][j]] / 3;
    nvsf.push_back(tricolor);
    if (i == 0) {
      vfmin = tricolor;
      vfmax = tricolor;
    } else if (tricolor < vfmin)
      vfmin = tricolor;
    else if (tricolor > vfmax)
      vfmax = tricolor;
  }

  vfmax -= vfmin;
  vector<float> nvsf2;
  for (unsigned int i = 0; i < nvsf.size(); i++)
    nvsf2.push_back((nvsf[i] - vfmin) / vfmax);

  return nvsf2;
}

float FEA2::ComputeStrainEnergy() {
  // sE = a' · K · a = a' · F
  vector<vector<float>> vvat;
  vector<float> vvati;
  for (unsigned int i = 0; i < vva.size(); i++) {
    vvati.push_back(vva[i][0]);
  }
  vvat.push_back(vvati);

  vector<vector<float>> vvsE = MultiplyMatricesEigen(vvat, vvf);
  sE = vvsE[0][0];

  if (sE < 0.0) sE = -sE;

  CurrentSE = sE;
  return sE;
}


float FEA2::NormalizeStrainEnergy() {
  int nEl = Ksize / 3;
  nsE = sE / nEl;

  emin_ = std::numeric_limits<float>::max();
  emax_ = std::numeric_limits<float>::min();
  for (unsigned int i = 0; i < ve_.size(); i++) {
    float en = ve_[i];
    if (en <= 0.0) continue;
    if (en > nsE) en = nsE;
    if (en < emin_) emin_ = en;
    if (en > emax_) emax_ = en;
  }

  return nsE;
}

void FEA2::UpdateForces() {
  unsigned int nadd = 3 * (vMPsXYZN_u.size() + vMPsXYZN_u2.size() -
                           vMPsXYZN_t.size() - vMPsXYZN_t2.size());
  for (unsigned int i = 0; i < nadd; i++) {
    vector<float> vfi = vector<float>(1, 0.0);
    vvf.push_back(vfi);
  }
}

void FEA2::ComputeNewDisplacement() {
  vva2 = MultiplyMatricesEigen(Ku, vvf);
  vva2 = vector_resize_cols(vva2, 3);
}

vector<vector<float>> FEA2::vector_resize_cols(vector<vector<float>> v1,
                                               unsigned int n) {
  vector<vector<float>> v2;
  vector<float> v2i;
  for (unsigned int i = 0; i < v1.size(); i++) {
    for (unsigned int j = 0; j < v1[i].size(); j++) {
      v2i.push_back(v1[i][j]);
      if (v2i.size() == n) {
        v2.push_back(v2i);
        v2i.clear();
      }
    }
  }
  return v2;
}

vector<vector<int>> FEA2::vector_resize_cols_int(vector<vector<int>> v1,
                                                 unsigned int n) {
  vector<vector<int>> v2;
  vector<int> v2i;
  for (unsigned int i = 0; i < v1.size(); i++) {
    for (unsigned int j = 0; j < v1[i].size(); j++) {
      v2i.push_back(v1[i][j]);
      if (v2i.size() == n) {
        v2.push_back(v2i);
        v2i.clear();
      }
    }
  }
  return v2;
}

float FEA2::GetStrainEnergy() { return sE; }

float FEA2::GetNormalizedStrainEnergy() { return nsE; }

void FEA2::ComputeInitialStrainEnergy(vector<vector<float>> vPoints) {
  Set_uf(vPoints);
  ComputeDisplacement();
  ComputeForces();
  sE0 = ComputeStrainEnergy();
  nsE0 = NormalizeStrainEnergy();
}

void FEA2::set_camera_pose(int model, Eigen::Vector4d qvec, Eigen::Vector3d tvec){
  qvecs_[model] = qvec;
  tvecs_[model] = tvec;
}

void FEA2::add_point3d(int model, Point3D *pMP) {
  if (model == 0 || model == 1) {
    vpMPs_t.push_back(pMP);
    vpMPs_pairs_[model].push_back(vpMPs_t.size() - 1);
  }
  else {
    vpMPs_ut.push_back(pMP);
  }
}

void FEA2::setbfea2(bool bSet) { bInFEA2 = bSet; }

void FEA2::setbfea(bool bSet) { bInFEA = bSet; }

void FEA2::DisplayPolygons(std::string &img_path,
                           std::vector<std::vector<int>> &tri,
                           std::vector<std::vector<float>> &points_2d,
                           std::vector<int> color) {
  cv::Mat img = cv::imread(img_path, cv::IMREAD_COLOR);

  std::vector<std::vector<cv::Point>> roi;

  for (unsigned int i = 0; i < tri.size(); i++) {
    cv::Point pt0 = cv::Point(points_2d[tri[i][0]][0], points_2d[tri[i][0]][1]);
    cv::Point pt1 = cv::Point(points_2d[tri[i][1]][0], points_2d[tri[i][1]][1]);
    cv::Point pt2 = cv::Point(points_2d[tri[i][2]][0], points_2d[tri[i][2]][1]);

    cv::line(img, pt0, pt1, cv::Scalar(color[0], color[1], color[2]), 2);
    cv::line(img, pt1, pt2, cv::Scalar(color[0], color[1], color[2]), 2);
    cv::line(img, pt2, pt0, cv::Scalar(color[0], color[1], color[2]), 2);

    std::vector<cv::Point> roi_i = {pt0, pt1, pt2};
    roi.push_back(roi_i);
  }

  img = SetTransparentColor(img, roi, 1.0, 0.5);

  cv::imshow("Image", img);

  cv::waitKey(0);

  cv::destroyAllWindows();
}

cv::Mat FEA2::SetTransparentColor(cv::Mat &img,
                                  std::vector<std::vector<cv::Point>> &roi,
                                  double alpha1, double alpha2) {
  cv::Mat out;
  cv::Mat layer = cv::Mat::zeros(img.size(), CV_8UC3);

  for (unsigned int i = 0; i < roi.size(); i++) {
    cv::Scalar colori = SetColor(vpMPs2DrawWgt[i]);

    std::vector<std::vector<cv::Point>> roi1;
    roi1.push_back(roi[i]);

    fillPoly(layer, roi1, colori);
  }

  cv::addWeighted(img, alpha1, layer, alpha2, 0, out);

  return out;
}

cv::Scalar FEA2::SetColor(float floatValue) {
  int cRed = floatValue * 255;
  int cBlue = 255 - floatValue * 255;
  int cGreen = 0;

  if (floatValue >= 0 && floatValue <= 0.5)
    cGreen = floatValue * 512;
  else if (floatValue > 0.5 && floatValue <= 1)
    cGreen = 255 - (floatValue - 0.5) * 512;
  else
    return -1;

  cv::Scalar output = cv::Scalar(cBlue, cGreen, cRed);
  return output;
}


void FEA2::set_mesh_settings(float mls_search_radius, int mls_polynomial_order,
                       float mesh_mu, float mesh_search_radius,
                       int mesh_max_neighbours, int mesh_surf_angle,
                       int mesh_min_angle,
                       int mesh_max_angle) {
  mls_search_radius_ = mls_search_radius;
  mls_polynomial_order_ = mls_polynomial_order;
  mesh_mu_ = mesh_mu;
  mesh_search_radius_ = mesh_search_radius;
  mesh_max_neighbours_ = mesh_max_neighbours;
  mesh_surf_angle_ = mesh_surf_angle;
  mesh_min_angle_ = mesh_min_angle;
  mesh_max_angle_ = mesh_max_angle;
}