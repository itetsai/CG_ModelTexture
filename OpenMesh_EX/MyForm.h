#pragma once

#include "DotNetUtilities.h"
#include "Mesh/GUA_OM.h"
#include "Mesh/DP.h"
#include <vector>
#include "math.h"
#include <Eigen/Sparse>

#include <iostream>
#include <stdlib.h>
#include <string>
#include <fstream>

#ifdef __APPLE__
#include <OpenGL/OpenGL.h>
#include <GLUT/glut.h>
#else
#include "../freeglut-3.0.0/include/glut.h"
#endif

#include "imageloader.h"

using namespace Eigen;

GLuint loadTexture(LImage* image) 
{
	GLuint textureId;
	glGenTextures(1, &textureId);
	glBindTexture(GL_TEXTURE_2D, textureId);
	glTexImage2D(GL_TEXTURE_2D,
		0,
		GL_RGB,
		image->width, image->height,
		0,
		GL_RGB,
		GL_UNSIGNED_BYTE,
		image->pixels);
	return textureId;
}

struct TexturePoint
{
	OMT::VertexHandle p;
	double x;
	double y;
};



struct EdgeWeight
{
	OMT::EdgeHandle edge;
	double weight;
};


#pragma region MESH
Tri_Mesh *mesh;
Tri_Mesh *mesh2;
std::vector<OMT::FVIter> fv;


int _textureId;
xform xf;
GLCamera camera;
float fov = 0.7f;

bool rrr = false;
int boundary_rotate = 0;
OMT::VIter pv;
std::vector<OMT::VFIter> f_it_list;
std::vector<TexturePoint> AllPoint;
#pragma endregion

#pragma region UI

struct FileData
{
	std::string filename;
	std::string filelocation;
};

struct Face
{
	std::vector<TexturePoint> vertex;
};

struct TextureMesh
{
	std::vector<TexturePoint> vertex;
	std::vector<Face> face;
	std::string textureName;
	std::string SetName;
	int textureID;
	int id;
};

std::vector<TextureMesh> All_texture_picture;
std::vector<TextureMesh> Show_texture;

std::vector<FileData> fd_list;

std::string picturename;

int setamount = 0;
#pragma endregion



double TwoPointLength(double x1, double y1, double z1, double x2, double y2, double z2)
{
	return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2));
}

double Holen(double a,double b,double c)
{
	double s = a + b + c;
	s = s*0.5;
	double qqo = s*(s - a)*(s - b)*(s - c);
	double nan = sqrt(qqo);
	return nan;
}


struct AssociationPoint
{
	OMT::VertexHandle p;
	std::vector<OMT::VertexHandle> association;
};

struct AssociationEdgeWeight
{
	OMT::VertexHandle p;
	std::vector<EdgeWeight> association;
};


static const Mouse::button physical_to_logical_map[] = {
	Mouse::NONE, Mouse::ROTATE, Mouse::MOVEXY, Mouse::MOVEZ,
	Mouse::MOVEZ, Mouse::MOVEXY, Mouse::MOVEXY, Mouse::MOVEXY,
};
Mouse::button Mouse_State = Mouse::ROTATE;


OMT::VertexHandle FindEdgePointToPoint(OMT::VertexHandle vh,OMT::EdgeHandle eh, Tri_Mesh m)
{
	using namespace OMT;

	HEHandle _hedge = m.halfedge_handle(eh, 0);
	VertexHandle v;

	v = m.from_vertex_handle(_hedge);
	if (v != vh)
	{
		return v;
	}

	v = m.to_vertex_handle(_hedge);
	if (v != vh)
	{
		return v;
	}
}

double ComputeWeight(OMT::EdgeHandle eh, Tri_Mesh m)
{
	using namespace OMT;

	HEHandle _hedge = m.halfedge_handle(eh, 0);
	VertexHandle v1 = m.from_vertex_handle(_hedge);
	VertexHandle v2 = m.to_vertex_handle(_hedge);
	vector<VertexHandle> two;
	for (VEIter ve = m.ve_iter(v1); ve; ve++)
	{
		VertexHandle fp1,fp2;
		fp1 = FindEdgePointToPoint(v1, ve.handle(), m);
		for (VEIter ve2 = m.ve_iter(fp1); ve2; ve2++)
		{
			fp2 = FindEdgePointToPoint(fp1, ve2.handle(), m);
			if (fp2 == v2)
				two.push_back(fp1);
		}
	}

	double cos1, cos2;

	double mid;
	mid = TwoPointLength(m.point(v1).data()[0], m.point(v1).data()[1], m.point(v1).data()[2], m.point(v2).data()[0], m.point(v2).data()[1], m.point(v2).data()[2]);

	double s1, s2;
	s1 = TwoPointLength(m.point(v1).data()[0], m.point(v1).data()[1], m.point(v1).data()[2], m.point(two[0]).data()[0], m.point(two[0]).data()[1], m.point(two[0]).data()[2]);
	s2 = TwoPointLength(m.point(v2).data()[0], m.point(v2).data()[1], m.point(v2).data()[2], m.point(two[0]).data()[0], m.point(two[0]).data()[1], m.point(two[0]).data()[2]);

	cos1 = (pow(s1, 2) + pow(s2, 2) - pow(mid, 2)) / (2 * s1*s2);

	s1 = TwoPointLength(m.point(v1).data()[0], m.point(v1).data()[1], m.point(v1).data()[2], m.point(two[1]).data()[0], m.point(two[1]).data()[1], m.point(two[1]).data()[2]);
	s2 = TwoPointLength(m.point(v2).data()[0], m.point(v2).data()[1], m.point(v2).data()[2], m.point(two[1]).data()[0], m.point(two[1]).data()[1], m.point(two[1]).data()[2]);

	cos2= (pow(s1, 2) + pow(s2, 2) - pow(mid, 2)) / (2 * s1*s2);
	
	double ct1 = 1 / tan(acos(cos2)), ct2 = 1 / tan(acos(cos1));

	double need = ct1 + ct2;
	if (need<0.00000001&&need>-0.00000001)
		need = 0;

	return need;
}


void findtotallength(std::vector<OMT::VertexHandle>&vh, OMT::VertexHandle now_vh, Tri_Mesh m, OMT::VertexHandle oringinal, OMT::VertexHandle last_vh)
{
	using namespace OMT;
	for (VEIter ve = m.ve_iter(now_vh); ve; ve++) //先找該點的edge
	{
		int face = 0;
		for (FIter fi = m.faces_begin(); fi != m.faces_end(); fi++)
		{
			for (FEIter fei = m.fe_iter(fi); fei; fei++)
			{
				if (ve.handle() == fei.handle())
					face++;//算與edge碰到的面有幾個
			}
		}

		if (face == 1)//與edge碰到的面只有一個代表是最外圍的edge
		{
			VertexHandle v = FindEdgePointToPoint(now_vh, ve.handle(), m);

			if (v != oringinal&& v != last_vh)
			{
				vh.push_back(v);
				findtotallength(vh, v, m, oringinal, now_vh);
				return;
			}
		}

	}
}
void firstboundary(std::vector<OMT::VertexHandle>&vh, OMT::VertexHandle now_vh, Tri_Mesh m, OMT::VertexHandle oringinal)
{
	using namespace OMT;
	for (VEIter ve = m.ve_iter(now_vh); ve; ve++) //先找該點的edge
	{
		int face = 0;
		for (FIter fi = m.faces_begin(); fi != m.faces_end(); fi++)
		{

			for (FEIter fei = m.fe_iter(fi); fei; fei++)
			{
				if (ve.handle() == fei.handle())
					face++;//算與edge碰到的面有幾個
			}
		}

		if (face == 1)//與edge碰到的面只有一個代表是最外圍的edge
		{
			VertexHandle v=FindEdgePointToPoint(now_vh,ve.handle(),m);

			if (v != oringinal)
			{
				vh.push_back(v);
				findtotallength(vh, v, m, oringinal, now_vh);
				return;
			}
		}

	}
}

void SortBoundary(Tri_Mesh m,std::vector<OMT::VertexHandle> &vh)
{
	using namespace OMT;
	
	double total_length;
	for (VIter vi = m.vertices_begin(); vi != m.vertices_end(); vi++)
	{
		if (m.is_boundary(vi))
		{			
			vh.push_back(vi);
			firstboundary(vh, vi, m,vi);
			break;
		}
	}
		
}

std::vector<TexturePoint> DefineSidePoint(std::vector<OMT::VertexHandle> &vh, Tri_Mesh m)//決定邊上的點
{
	double total_length = 0;
	for (int i = 0; i < vh.size()-1; i++)
	{
		total_length += TwoPointLength(m.point(vh[i]).data()[0], m.point(vh[i]).data()[1], m.point(vh[i]).data()[2], m.point(vh[i + 1]).data()[0], m.point(vh[i + 1]).data()[1], m.point(vh[i + 1]).data()[2]);
	}
	total_length += TwoPointLength(m.point(vh[vh.size() - 1]).data()[0], m.point(vh[vh.size() - 1]).data()[1], m.point(vh[vh.size() - 1]).data()[2], m.point(vh[0]).data()[0], m.point(vh[0]).data()[1], m.point(vh[0]).data()[2]);
	total_length /= 4.0;

	std::vector<TexturePoint>tp(vh.size());
	double now_length = 0;

	for (int i = 0; i < vh.size(); i++)
	{
		tp[i].p = vh[i];
	}

	tp[0].x = 0; tp[0].y = 0;

	for (int i = 0; i < vh.size()-1; i++)
	{
		now_length += TwoPointLength(m.point(vh[i]).data()[0], m.point(vh[i]).data()[1], m.point(vh[i]).data()[2], m.point(vh[i + 1]).data()[0], m.point(vh[i + 1]).data()[1], m.point(vh[i + 1]).data()[2])/total_length;
		if (now_length < 1)//左
		{
			tp[i + 1].x = 0;
			tp[i + 1].y = now_length;
		}
		else if (now_length >= 1 && now_length < 2)//上
		{
			tp[i + 1].x = now_length - 1.0;
			tp[i + 1].y = 1;
		}
		else if (now_length >= 2 && now_length < 3)//右
		{
			tp[i + 1].x = 1;
			tp[i + 1].y = 1 - (now_length - 2.0);
		}
		else if (now_length >= 3 && now_length < 4)//下
		{
			tp[i + 1].x = 1 - (now_length - 3.0);
			tp[i + 1].y = 0;
		}
	}

	return tp;
}

std::vector<TexturePoint> Advanced_Fill_Matrix(Tri_Mesh m, std::vector<TexturePoint> boundary_list, std::vector<EdgeWeight> weight_list)
{
	using namespace OMT;
	int size = 0;


	vector<VertexHandle> vh;//內部點的矩陣
	for (OMT::VIter v_it = m.vertices_begin(); v_it != m.vertices_end(); ++v_it)//算有幾個內部點
	{
		if (!m.is_boundary(v_it.handle()))
		{
			size++;
			vh.push_back(v_it.handle());//將內部點加入矩陣
		}
	}


	if (size == 0)
	{
		vector<TexturePoint> zero;
		return zero;
	}

	SparseMatrix<double> A(size, size);//決定矩陣大小
	vector<AssociationEdgeWeight> ae_list;

	for (int i = 0; i < vh.size(); i++)//建立每個點相關的edge
	{
		AssociationEdgeWeight ae;
		ae.p = vh[i];
		for (VEIter ve = m.ve_iter(vh[i]); ve; ve++)
		{
			EdgeWeight ew;

			ew.edge = ve.handle();
			for (int j = 0; j < weight_list.size(); j++)
				if (ew.edge == weight_list[j].edge)
					ew.weight = weight_list[j].weight;
			ae.association.push_back(ew);		
		}
		ae_list.push_back(ae);
	}

	for (int i = 0; i < ae_list.size(); i++)//算內部點的權重
	{
		double w = 0;
		for (int j = 0; j < ae_list[i].association.size(); j++)
			w += ae_list[i].association[j].weight;

		A.insert(i, i) = w;
	}


	for (int i = 0; i<vh.size(); i++)
		for (int j = 0; j < vh.size(); j++)
			if (i != j)
			{
				double w=0;
				for (int k = 0; k < ae_list[i].association.size(); k++)
					for (int m = 0; m < ae_list[j].association.size(); m++)
						if (ae_list[i].association[k].edge == ae_list[j].association[m].edge)
							w = ae_list[i].association[k].weight;
				
				A.insert(i, j) = -1*w;
			}
			
	A.makeCompressed();
//A矩陣填完/////////////////////////////////////////////////////////////////////////////////


	std::vector<VectorXd> B;
	B.resize(2);

	B[0].resize(vh.size());
	B[1].resize(vh.size());

	for (int i = 0; i < vh.size(); i++)
	{
		double wx = 0;
		double wy = 0;
		for (int j = 0; j < ae_list[i].association.size(); j++)
		{
			VertexHandle v;
			v = FindEdgePointToPoint(ae_list[i].p, ae_list[i].association[j].edge,m);
			if (m.is_boundary(v))
			{
				double w= ae_list[i].association[j].weight;
				TexturePoint tp;

				for (int k = 0; k < boundary_list.size(); k++)
				{
					if (v == boundary_list[k].p)
					{
						wx += w*boundary_list[k].x;
						wy += w*boundary_list[k].y;
					}
				}
			}
		}
		B[0][i] = wx;
		B[1][i] = wy;
	}
	//B矩陣填完////////////////////////////////////////////////////////////////////////////////////////////

	// 利用 SparseQR 來解 AX = B
	SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> linearSolver;
	linearSolver.compute(A);

	std::vector<VectorXd> X;
	X.resize(2);

	X[0] = linearSolver.solve(B[0]);
	X[1] = linearSolver.solve(B[1]);

	vector<TexturePoint> tp_list;

	for (int i = 0; i < vh.size(); i++)
	{
		TexturePoint tp;
		tp.p = ae_list[i].p;
		tp.x = X[0][i];
		tp.y = X[1][i];
		tp_list.push_back(tp);
	}//完成///////////////////////////////////////////////////
	return tp_list;
}


void SaveSetFile(std::string _fileName, std::vector<TextureMesh> save)
{
	using namespace std;

	string filename=_fileName+".txt";
	fstream file;
	file.open(filename, ios::out);
	
	file << save.size() << endl;
	for (int i = 0; i < save.size(); i++)
	{
		file << save[i].vertex.size() << endl;
		for (int j = 0; j < save[i].vertex.size(); j++)
			file << save[i].vertex[j].p << " " << save[i].vertex[j].x << " " << save[i].vertex[j].y << endl;
		
		file << endl;
		file << save[i].face.size() << endl;//

		for (int j = 0; j < save[i].face.size(); j++)
		{
			file << save[i].face[j].vertex[0].p << " " << save[i].face[j].vertex[0].x << " " << save[i].face[j].vertex[0].y << endl;
			file << save[i].face[j].vertex[1].p << " " << save[i].face[j].vertex[1].x << " " << save[i].face[j].vertex[1].y << endl;
			file << save[i].face[j].vertex[2].p << " " << save[i].face[j].vertex[2].x << " " << save[i].face[j].vertex[2].y << endl;
			file << endl;
		}
		file << endl;
		file << save[i].textureName << endl;
		file << save[i].SetName << endl;
		file << save[i].id << endl;
	}

	file.close();
}

void LoadSetFile(std::string _fileName, std::vector<TextureMesh>& load)
{
	using namespace std;
	fstream file;
	file.open(_fileName, ios::in);
	int size;
	file >> size;

	for (int i = 0; i < size; i++)
	{
		int vertexsize;
		TextureMesh tm;
		file >> vertexsize;
		for (int j = 0; j < vertexsize; j++)
		{
			TexturePoint tp;
			int id;
			file >> id;
			tp.p.invalidate();
			tp.p.__increment();
			tp.p.__increment(id);
			file >> tp.x;
			file >> tp.y;
			tm.vertex.push_back(tp);
		}
		int facesize;
		file >> facesize;

		for (int j = 0; j < facesize; j++)
		{
			Face f;
			for (int k = 0; k < 3; k++)
			{
				TexturePoint tp;
				int id;
				file >> id;
				tp.p.invalidate();
				tp.p.__increment();
				tp.p.__increment(id);
				file >> tp.x;
				file >> tp.y;
				f.vertex.push_back(tp);
			}
			tm.face.push_back(f);
		}

		file >> tm.textureName;
		file >> tm.SetName;
		file >> tm.id;
		load.push_back(tm);
	}
	file.close();
}

namespace OpenMesh_EX {

	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;
	using namespace OMT;
	/// <summary>
	/// MyForm 的摘要
	/// </summary>
	public ref class MyForm : public System::Windows::Forms::Form
	{
#pragma region dontchange

	public:
		MyForm(void)
		{
			InitializeComponent();
			//
			//TODO:  在此加入建構函式程式碼
			//
		}

	protected:
		/// <summary>
		/// 清除任何使用中的資源。
		/// </summary>
		~MyForm()
		{
			if (components)
			{
				delete components;
			}
		}

	private: System::Windows::Forms::MenuStrip^  menuStrip1;
	private: System::Windows::Forms::ToolStripMenuItem^  fileToolStripMenuItem;
	private: System::Windows::Forms::ToolStripMenuItem^  loadModelToolStripMenuItem;
	private: System::Windows::Forms::OpenFileDialog^  openModelDialog;
	private: System::Windows::Forms::SaveFileDialog^  saveModelDialog;
	private: System::Windows::Forms::ToolStripMenuItem^  saveModelToolStripMenuItem;
	private: HKOGLPanel::HKOGLPanelControl^  hkoglPanelControl1;
	private: System::Windows::Forms::RadioButton^  radioButton1;
	private: System::Windows::Forms::RadioButton^  radioButton2;
	private: System::Windows::Forms::Button^  button1;
	private: System::Windows::Forms::Button^  button2;
	private: System::Windows::Forms::Button^  button3;
	private: System::Windows::Forms::PictureBox^  pictureBox1;
	private: System::Windows::Forms::ListBox^  listBox1;
	private: System::Windows::Forms::Panel^  panel1;
	private: System::Windows::Forms::Panel^  panel2;
	private: System::Windows::Forms::RadioButton^  radioButton4;
	private: System::Windows::Forms::RadioButton^  radioButton3;
	private: System::Windows::Forms::OpenFileDialog^  PictureLoad;
	private: System::Windows::Forms::Button^  button4;

	private: System::Windows::Forms::SaveFileDialog^  SaveSet;
	private: System::Windows::Forms::Button^  button5;
	private: System::Windows::Forms::ListBox^  listBox2;
	private: System::Windows::Forms::CheckBox^  checkBox1;
	private: System::Windows::Forms::TextBox^  textBox1;
	private: System::Windows::Forms::Button^  button6;
	private: System::Windows::Forms::OpenFileDialog^  LoadSet;
	private: System::Windows::Forms::Label^  label1;
	private: System::Windows::Forms::Label^  label2;
	private: System::Windows::Forms::Button^  button7;

	protected:

	private:
		/// <summary>
		/// 設計工具所需的變數。
		/// </summary>
		System::ComponentModel::Container ^components;

#pragma region Windows Form Designer generated code
		/// <summary>
		/// 此為設計工具支援所需的方法 - 請勿使用程式碼編輯器修改
		/// 這個方法的內容。
		/// </summary>
		void InitializeComponent(void)
		{
			HKOGLPanel::HKCOGLPanelCameraSetting^  hkcoglPanelCameraSetting1 = (gcnew HKOGLPanel::HKCOGLPanelCameraSetting());
			HKOGLPanel::HKCOGLPanelPixelFormat^  hkcoglPanelPixelFormat1 = (gcnew HKOGLPanel::HKCOGLPanelPixelFormat());
			this->menuStrip1 = (gcnew System::Windows::Forms::MenuStrip());
			this->fileToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->loadModelToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->saveModelToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->openModelDialog = (gcnew System::Windows::Forms::OpenFileDialog());
			this->saveModelDialog = (gcnew System::Windows::Forms::SaveFileDialog());
			this->hkoglPanelControl1 = (gcnew HKOGLPanel::HKOGLPanelControl());
			this->radioButton1 = (gcnew System::Windows::Forms::RadioButton());
			this->radioButton2 = (gcnew System::Windows::Forms::RadioButton());
			this->button1 = (gcnew System::Windows::Forms::Button());
			this->button2 = (gcnew System::Windows::Forms::Button());
			this->button3 = (gcnew System::Windows::Forms::Button());
			this->pictureBox1 = (gcnew System::Windows::Forms::PictureBox());
			this->listBox1 = (gcnew System::Windows::Forms::ListBox());
			this->panel1 = (gcnew System::Windows::Forms::Panel());
			this->panel2 = (gcnew System::Windows::Forms::Panel());
			this->radioButton4 = (gcnew System::Windows::Forms::RadioButton());
			this->radioButton3 = (gcnew System::Windows::Forms::RadioButton());
			this->PictureLoad = (gcnew System::Windows::Forms::OpenFileDialog());
			this->button4 = (gcnew System::Windows::Forms::Button());
			this->SaveSet = (gcnew System::Windows::Forms::SaveFileDialog());
			this->button5 = (gcnew System::Windows::Forms::Button());
			this->listBox2 = (gcnew System::Windows::Forms::ListBox());
			this->checkBox1 = (gcnew System::Windows::Forms::CheckBox());
			this->textBox1 = (gcnew System::Windows::Forms::TextBox());
			this->button6 = (gcnew System::Windows::Forms::Button());
			this->LoadSet = (gcnew System::Windows::Forms::OpenFileDialog());
			this->label1 = (gcnew System::Windows::Forms::Label());
			this->label2 = (gcnew System::Windows::Forms::Label());
			this->button7 = (gcnew System::Windows::Forms::Button());
			this->menuStrip1->SuspendLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pictureBox1))->BeginInit();
			this->panel1->SuspendLayout();
			this->panel2->SuspendLayout();
			this->SuspendLayout();
			// 
			// menuStrip1
			// 
			this->menuStrip1->ImageScalingSize = System::Drawing::Size(20, 20);
			this->menuStrip1->Items->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(1) { this->fileToolStripMenuItem });
			this->menuStrip1->Location = System::Drawing::Point(0, 0);
			this->menuStrip1->Name = L"menuStrip1";
			this->menuStrip1->Padding = System::Windows::Forms::Padding(8, 2, 0, 2);
			this->menuStrip1->Size = System::Drawing::Size(1450, 27);
			this->menuStrip1->TabIndex = 1;
			this->menuStrip1->Text = L"menuStrip1";
			// 
			// fileToolStripMenuItem
			// 
			this->fileToolStripMenuItem->DropDownItems->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(2) {
				this->loadModelToolStripMenuItem,
					this->saveModelToolStripMenuItem
			});
			this->fileToolStripMenuItem->Name = L"fileToolStripMenuItem";
			this->fileToolStripMenuItem->Size = System::Drawing::Size(45, 23);
			this->fileToolStripMenuItem->Text = L"File";
			// 
			// loadModelToolStripMenuItem
			// 
			this->loadModelToolStripMenuItem->Name = L"loadModelToolStripMenuItem";
			this->loadModelToolStripMenuItem->Size = System::Drawing::Size(168, 26);
			this->loadModelToolStripMenuItem->Text = L"Load Model";
			this->loadModelToolStripMenuItem->Click += gcnew System::EventHandler(this, &MyForm::loadModelToolStripMenuItem_Click);
			// 
			// saveModelToolStripMenuItem
			// 
			this->saveModelToolStripMenuItem->Name = L"saveModelToolStripMenuItem";
			this->saveModelToolStripMenuItem->Size = System::Drawing::Size(168, 26);
			this->saveModelToolStripMenuItem->Text = L"Save Model";
			this->saveModelToolStripMenuItem->Click += gcnew System::EventHandler(this, &MyForm::saveModelToolStripMenuItem_Click);
			// 
			// openModelDialog
			// 
			this->openModelDialog->FileOk += gcnew System::ComponentModel::CancelEventHandler(this, &MyForm::openModelDialog_FileOk);
			// 
			// saveModelDialog
			// 
			this->saveModelDialog->DefaultExt = L"obj";
			this->saveModelDialog->FileOk += gcnew System::ComponentModel::CancelEventHandler(this, &MyForm::saveModelDialog_FileOk);
			// 
			// hkoglPanelControl1
			// 
			this->hkoglPanelControl1->Anchor = static_cast<System::Windows::Forms::AnchorStyles>(((System::Windows::Forms::AnchorStyles::Top | System::Windows::Forms::AnchorStyles::Bottom)
				| System::Windows::Forms::AnchorStyles::Left));
			hkcoglPanelCameraSetting1->Far = 1000;
			hkcoglPanelCameraSetting1->Fov = 45;
			hkcoglPanelCameraSetting1->Near = -1000;
			hkcoglPanelCameraSetting1->Type = HKOGLPanel::HKCOGLPanelCameraSetting::CAMERATYPE::ORTHOGRAPHIC;
			this->hkoglPanelControl1->Camera_Setting = hkcoglPanelCameraSetting1;
			this->hkoglPanelControl1->Location = System::Drawing::Point(0, 31);
			this->hkoglPanelControl1->Margin = System::Windows::Forms::Padding(4);
			this->hkoglPanelControl1->Name = L"hkoglPanelControl1";
			hkcoglPanelPixelFormat1->Accumu_Buffer_Bits = HKOGLPanel::HKCOGLPanelPixelFormat::PIXELBITS::BITS_0;
			hkcoglPanelPixelFormat1->Alpha_Buffer_Bits = HKOGLPanel::HKCOGLPanelPixelFormat::PIXELBITS::BITS_0;
			hkcoglPanelPixelFormat1->Stencil_Buffer_Bits = HKOGLPanel::HKCOGLPanelPixelFormat::PIXELBITS::BITS_0;
			this->hkoglPanelControl1->Pixel_Format = hkcoglPanelPixelFormat1;
			this->hkoglPanelControl1->Size = System::Drawing::Size(817, 578);
			this->hkoglPanelControl1->TabIndex = 2;
			this->hkoglPanelControl1->Load += gcnew System::EventHandler(this, &MyForm::hkoglPanelControl1_Load);
			this->hkoglPanelControl1->Paint += gcnew System::Windows::Forms::PaintEventHandler(this, &MyForm::hkoglPanelControl1_Paint);
			this->hkoglPanelControl1->MouseDown += gcnew System::Windows::Forms::MouseEventHandler(this, &MyForm::hkoglPanelControl1_MouseDown);
			this->hkoglPanelControl1->MouseMove += gcnew System::Windows::Forms::MouseEventHandler(this, &MyForm::hkoglPanelControl1_MouseMove);
			this->hkoglPanelControl1->MouseWheel += gcnew System::Windows::Forms::MouseEventHandler(this, &MyForm::hkoglPanelControl1_MouseWheel);
			// 
			// radioButton1
			// 
			this->radioButton1->AutoSize = true;
			this->radioButton1->Location = System::Drawing::Point(3, 16);
			this->radioButton1->Margin = System::Windows::Forms::Padding(3, 2, 3, 2);
			this->radioButton1->Name = L"radioButton1";
			this->radioButton1->Size = System::Drawing::Size(58, 19);
			this->radioButton1->TabIndex = 3;
			this->radioButton1->Text = L"選面";
			this->radioButton1->UseVisualStyleBackColor = true;
			this->radioButton1->CheckedChanged += gcnew System::EventHandler(this, &MyForm::radioButton1_CheckedChanged);
			// 
			// radioButton2
			// 
			this->radioButton2->AutoSize = true;
			this->radioButton2->Checked = true;
			this->radioButton2->Location = System::Drawing::Point(3, 39);
			this->radioButton2->Margin = System::Windows::Forms::Padding(3, 2, 3, 2);
			this->radioButton2->Name = L"radioButton2";
			this->radioButton2->Size = System::Drawing::Size(58, 19);
			this->radioButton2->TabIndex = 4;
			this->radioButton2->TabStop = true;
			this->radioButton2->Text = L"旋轉";
			this->radioButton2->UseVisualStyleBackColor = true;
			this->radioButton2->CheckedChanged += gcnew System::EventHandler(this, &MyForm::radioButton2_CheckedChanged);
			// 
			// button1
			// 
			this->button1->Location = System::Drawing::Point(851, 165);
			this->button1->Margin = System::Windows::Forms::Padding(3, 2, 3, 2);
			this->button1->Name = L"button1";
			this->button1->Size = System::Drawing::Size(127, 29);
			this->button1->TabIndex = 5;
			this->button1->Text = L"Mapping";
			this->button1->UseVisualStyleBackColor = true;
			this->button1->Click += gcnew System::EventHandler(this, &MyForm::button1_Click);
			// 
			// button2
			// 
			this->button2->Location = System::Drawing::Point(852, 200);
			this->button2->Margin = System::Windows::Forms::Padding(4);
			this->button2->Name = L"button2";
			this->button2->Size = System::Drawing::Size(127, 29);
			this->button2->TabIndex = 6;
			this->button2->Text = L"Rotate Boundary";
			this->button2->UseVisualStyleBackColor = true;
			this->button2->Click += gcnew System::EventHandler(this, &MyForm::button2_Click);
			// 
			// button3
			// 
			this->button3->Location = System::Drawing::Point(877, 549);
			this->button3->Name = L"button3";
			this->button3->Size = System::Drawing::Size(96, 33);
			this->button3->TabIndex = 7;
			this->button3->Text = L"Load Picture";
			this->button3->UseVisualStyleBackColor = true;
			this->button3->Click += gcnew System::EventHandler(this, &MyForm::button3_Click);
			// 
			// pictureBox1
			// 
			this->pictureBox1->BackColor = System::Drawing::SystemColors::ControlLight;
			this->pictureBox1->ImageLocation = L"";
			this->pictureBox1->Location = System::Drawing::Point(1015, 47);
			this->pictureBox1->Name = L"pictureBox1";
			this->pictureBox1->Size = System::Drawing::Size(316, 300);
			this->pictureBox1->SizeMode = System::Windows::Forms::PictureBoxSizeMode::Zoom;
			this->pictureBox1->TabIndex = 8;
			this->pictureBox1->TabStop = false;
			// 
			// listBox1
			// 
			this->listBox1->FormattingEnabled = true;
			this->listBox1->ItemHeight = 15;
			this->listBox1->Location = System::Drawing::Point(1176, 397);
			this->listBox1->Name = L"listBox1";
			this->listBox1->Size = System::Drawing::Size(155, 169);
			this->listBox1->TabIndex = 9;
			this->listBox1->SelectedIndexChanged += gcnew System::EventHandler(this, &MyForm::listBox1_SelectedIndexChanged);
			// 
			// panel1
			// 
			this->panel1->Controls->Add(this->radioButton2);
			this->panel1->Controls->Add(this->radioButton1);
			this->panel1->Location = System::Drawing::Point(851, 97);
			this->panel1->Name = L"panel1";
			this->panel1->Size = System::Drawing::Size(97, 63);
			this->panel1->TabIndex = 10;
			// 
			// panel2
			// 
			this->panel2->Controls->Add(this->radioButton4);
			this->panel2->Controls->Add(this->radioButton3);
			this->panel2->Location = System::Drawing::Point(851, 236);
			this->panel2->Name = L"panel2";
			this->panel2->Size = System::Drawing::Size(122, 72);
			this->panel2->TabIndex = 11;
			// 
			// radioButton4
			// 
			this->radioButton4->AutoSize = true;
			this->radioButton4->Checked = true;
			this->radioButton4->Location = System::Drawing::Point(3, 35);
			this->radioButton4->Name = L"radioButton4";
			this->radioButton4->Size = System::Drawing::Size(89, 19);
			this->radioButton4->TabIndex = 1;
			this->radioButton4->TabStop = true;
			this->radioButton4->Text = L"Wireframe";
			this->radioButton4->UseVisualStyleBackColor = true;
			this->radioButton4->CheckedChanged += gcnew System::EventHandler(this, &MyForm::radioButton4_CheckedChanged);
			// 
			// radioButton3
			// 
			this->radioButton3->AutoSize = true;
			this->radioButton3->Location = System::Drawing::Point(3, 10);
			this->radioButton3->Name = L"radioButton3";
			this->radioButton3->Size = System::Drawing::Size(72, 19);
			this->radioButton3->TabIndex = 0;
			this->radioButton3->Text = L"Texture";
			this->radioButton3->UseVisualStyleBackColor = true;
			this->radioButton3->CheckedChanged += gcnew System::EventHandler(this, &MyForm::radioButton3_CheckedChanged);
			// 
			// PictureLoad
			// 
			this->PictureLoad->FileName = L"openFileDialog1";
			this->PictureLoad->FileOk += gcnew System::ComponentModel::CancelEventHandler(this, &MyForm::PictureLoad_FileOk);
			// 
			// button4
			// 
			this->button4->Location = System::Drawing::Point(932, 390);
			this->button4->Name = L"button4";
			this->button4->Size = System::Drawing::Size(75, 25);
			this->button4->TabIndex = 12;
			this->button4->Text = L"Setting";
			this->button4->UseVisualStyleBackColor = true;
			this->button4->Click += gcnew System::EventHandler(this, &MyForm::button4_Click);
			// 
			// SaveSet
			// 
			this->SaveSet->FileOk += gcnew System::ComponentModel::CancelEventHandler(this, &MyForm::SaveSet_FileOk);
			// 
			// button5
			// 
			this->button5->Location = System::Drawing::Point(877, 510);
			this->button5->Name = L"button5";
			this->button5->Size = System::Drawing::Size(96, 33);
			this->button5->TabIndex = 13;
			this->button5->Text = L"Save Set";
			this->button5->UseVisualStyleBackColor = true;
			this->button5->Click += gcnew System::EventHandler(this, &MyForm::button5_Click);
			// 
			// listBox2
			// 
			this->listBox2->FormattingEnabled = true;
			this->listBox2->ItemHeight = 15;
			this->listBox2->Location = System::Drawing::Point(1016, 397);
			this->listBox2->Name = L"listBox2";
			this->listBox2->Size = System::Drawing::Size(155, 169);
			this->listBox2->TabIndex = 14;
			this->listBox2->SelectedIndexChanged += gcnew System::EventHandler(this, &MyForm::listBox2_SelectedIndexChanged);
			// 
			// checkBox1
			// 
			this->checkBox1->AutoSize = true;
			this->checkBox1->Location = System::Drawing::Point(1015, 375);
			this->checkBox1->Name = L"checkBox1";
			this->checkBox1->Size = System::Drawing::Size(77, 19);
			this->checkBox1->TabIndex = 15;
			this->checkBox1->Text = L"complex";
			this->checkBox1->UseVisualStyleBackColor = true;
			this->checkBox1->CheckedChanged += gcnew System::EventHandler(this, &MyForm::checkBox1_CheckedChanged);
			// 
			// textBox1
			// 
			this->textBox1->Location = System::Drawing::Point(826, 390);
			this->textBox1->Name = L"textBox1";
			this->textBox1->Size = System::Drawing::Size(100, 25);
			this->textBox1->TabIndex = 16;
			// 
			// button6
			// 
			this->button6->Location = System::Drawing::Point(877, 471);
			this->button6->Name = L"button6";
			this->button6->Size = System::Drawing::Size(96, 33);
			this->button6->TabIndex = 17;
			this->button6->Text = L"Load Set";
			this->button6->UseVisualStyleBackColor = true;
			this->button6->Click += gcnew System::EventHandler(this, &MyForm::button6_Click);
			// 
			// LoadSet
			// 
			this->LoadSet->FileName = L"openFileDialog1";
			this->LoadSet->FileOk += gcnew System::ComponentModel::CancelEventHandler(this, &MyForm::LoadSet_FileOk);
			// 
			// label1
			// 
			this->label1->AutoSize = true;
			this->label1->Font = (gcnew System::Drawing::Font(L"新細明體", 9, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(136)));
			this->label1->Location = System::Drawing::Point(1146, 379);
			this->label1->Name = L"label1";
			this->label1->Size = System::Drawing::Size(25, 15);
			this->label1->TabIndex = 18;
			this->label1->Text = L"Set";
			// 
			// label2
			// 
			this->label2->AutoSize = true;
			this->label2->Font = (gcnew System::Drawing::Font(L"新細明體", 9, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(136)));
			this->label2->Location = System::Drawing::Point(1284, 379);
			this->label2->Name = L"label2";
			this->label2->Size = System::Drawing::Size(47, 15);
			this->label2->TabIndex = 19;
			this->label2->Text = L"Picture";
			// 
			// button7
			// 
			this->button7->Location = System::Drawing::Point(877, 432);
			this->button7->Name = L"button7";
			this->button7->Size = System::Drawing::Size(96, 33);
			this->button7->TabIndex = 20;
			this->button7->Text = L"All Mapping";
			this->button7->UseVisualStyleBackColor = true;
			this->button7->Click += gcnew System::EventHandler(this, &MyForm::button7_Click);
			// 
			// MyForm
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(8, 15);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->ClientSize = System::Drawing::Size(1450, 609);
			this->Controls->Add(this->button7);
			this->Controls->Add(this->label2);
			this->Controls->Add(this->label1);
			this->Controls->Add(this->button6);
			this->Controls->Add(this->textBox1);
			this->Controls->Add(this->checkBox1);
			this->Controls->Add(this->listBox2);
			this->Controls->Add(this->button5);
			this->Controls->Add(this->button4);
			this->Controls->Add(this->panel2);
			this->Controls->Add(this->panel1);
			this->Controls->Add(this->listBox1);
			this->Controls->Add(this->pictureBox1);
			this->Controls->Add(this->button3);
			this->Controls->Add(this->button2);
			this->Controls->Add(this->button1);
			this->Controls->Add(this->hkoglPanelControl1);
			this->Controls->Add(this->menuStrip1);
			this->MainMenuStrip = this->menuStrip1;
			this->Margin = System::Windows::Forms::Padding(4);
			this->Name = L"MyForm";
			this->Text = L"OpenMesh_EX";
			this->Load += gcnew System::EventHandler(this, &MyForm::MyForm_Load);
			this->menuStrip1->ResumeLayout(false);
			this->menuStrip1->PerformLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pictureBox1))->EndInit();
			this->panel1->ResumeLayout(false);
			this->panel1->PerformLayout();
			this->panel2->ResumeLayout(false);
			this->panel2->PerformLayout();
			this->ResumeLayout(false);
			this->PerformLayout();

		}
#pragma endregion
#pragma endregion

private: System::Void hkoglPanelControl1_Load(System::Object^  sender, System::EventArgs^  e)
{

}


private: System::Void hkoglPanelControl1_Paint(System::Object^  sender, System::Windows::Forms::PaintEventArgs^  e)
{
	glEnable(GL_COLOR_MATERIAL);
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT);

	point center;
	center[0] = 0.0;
	center[1] = 0.0;
	center[2] = 0.0;
	camera.setupGL(xf * center, 1.0);

	glPushMatrix();
	glMatrixMode(GL_MODELVIEW);
	glMultMatrixd((double *)xf);

	for (int i = 0; i < Show_texture.size(); i++)///已貼的貼圖
	{
		glPushMatrix();

		glEnable(GL_LIGHTING);
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, Show_texture[i].textureID);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glBegin(GL_TRIANGLES);
		glColor3f(1.0, 1.0, 1.0);
		for (int j = 0; j < Show_texture[i].face.size(); j++)
		{
			glTexCoord2f(Show_texture[i].face[j].vertex[0].x, Show_texture[i].face[j].vertex[0].y);
			glNormal3dv(mesh->point(Show_texture[i].face[j].vertex[0].p).data());//不打光圖形有些地方會暗暗的
			glVertex3dv(mesh->point(Show_texture[i].face[j].vertex[0].p).data());

			glTexCoord2f(Show_texture[i].face[j].vertex[1].x, Show_texture[i].face[j].vertex[1].y);
			glNormal3dv(mesh->point(Show_texture[i].face[j].vertex[1].p).data());//不打光圖形有些地方會暗暗的
			glVertex3dv(mesh->point(Show_texture[i].face[j].vertex[1].p).data());

			glTexCoord2f(Show_texture[i].face[j].vertex[2].x, Show_texture[i].face[j].vertex[2].y);
			glNormal3dv(mesh->point(Show_texture[i].face[j].vertex[2].p).data());//不打光圖形有些地方會暗暗的
			glVertex3dv(mesh->point(Show_texture[i].face[j].vertex[2].p).data());
		}
		glEnd();
		glDisable(GL_TEXTURE_2D);
		glPopMatrix();
	}
	if (radioButton4->Checked)
	{
		for (int i = 0; i < Show_texture.size(); i++)///已貼的貼圖
		{
			glPushMatrix();
			glPointSize(8.0);
			glBegin(GL_POINTS);
			glColor3f(1, 0, 0);

			for (int j = 0; j < Show_texture[i].vertex.size(); j++)
				glVertex3dv(mesh->point(Show_texture[i].vertex[j].p).data());

			glEnd();
			glPopMatrix();
		}

		for (int i = 0; i < Show_texture.size(); i++)///已貼的貼圖
		{
			glPushMatrix();
			glDisable(GL_LIGHTING);
			glLineWidth(2.0);
			glBegin(GL_LINES);
			glColor3f(1.0, 0.3, 0);
			for (int j = 0; j < Show_texture[i].face.size(); j++)
			{

				glVertex3dv(mesh->point(Show_texture[i].face[j].vertex[0].p).data());
				glVertex3dv(mesh->point(Show_texture[i].face[j].vertex[1].p).data());

				glVertex3dv(mesh->point(Show_texture[i].face[j].vertex[1].p).data());
				glVertex3dv(mesh->point(Show_texture[i].face[j].vertex[2].p).data());

				glVertex3dv(mesh->point(Show_texture[i].face[j].vertex[2].p).data());
				glVertex3dv(mesh->point(Show_texture[i].face[j].vertex[0].p).data());
			}
			glEnd();
			glPopMatrix();
		}
	}

	if (rrr == true && mesh2 != NULL)
	{


		glPushMatrix();

		glEnable(GL_LIGHTING);



		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, _textureId);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glBegin(GL_TRIANGLES);
		glColor3f(1.0, 1.0, 1.0);
		for (OMT::FIter f_it = mesh2->faces_begin(); f_it != mesh2->faces_end(); ++f_it)
		{
			for (OMT::FVIter fv_it = mesh2->fv_iter(f_it); fv_it; ++fv_it)
			{
				double x, y;
				for (int i = 0; i < AllPoint.size(); i++)
					if (fv_it.handle() == AllPoint[i].p)
					{
						x = AllPoint[i].x;
						y = AllPoint[i].y;
					}
				glTexCoord2f(x, y);
				glNormal3dv(mesh2->point(fv_it.handle()).data());//不打光圖形有些地方會暗暗的
				glVertex3dv(mesh2->point(fv_it.handle()).data());

			}
		}

		glEnd();
		glDisable(GL_TEXTURE_2D);
		glPopMatrix();

	}

	if (mesh != NULL)
	{
		glPushMatrix();
		glEnable(GL_LIGHTING);
		glBegin(GL_TRIANGLES);
		glColor3f(0.0, 1.0, 0.0);
		for (int i = 0; i < f_it_list.size(); i++)
		{
			for (OMT::FVIter fv_it = mesh->fv_iter(f_it_list[i]); fv_it; ++fv_it)
			{
				glNormal3dv(mesh->point(fv_it.handle()).data());
				glVertex3dv(mesh->point(fv_it.handle()).data());
			}
		}

		glEnd();
		glPopMatrix();


		glPushMatrix();
		if (radioButton3->Checked)
			mesh->Render_Solid();

		if (radioButton4->Checked)
			mesh->Render_SolidWireframe();

		glPopMatrix();
	}


	glPopMatrix();
}

private: System::Void hkoglPanelControl1_MouseDown(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e)
{
	if (e->Button == System::Windows::Forms::MouseButtons::Left ||
		e->Button == System::Windows::Forms::MouseButtons::Middle)
	{
		point center;
		Mouse_State = Mouse::NONE;
		center[0] = 0.0;
		center[1] = 0.0;
		center[2] = 0.0;
		camera.mouse(e->X, e->Y, Mouse_State,
			xf * center,
			1.0, xf);

	}


	if (e->Button == System::Windows::Forms::MouseButtons::Left&&radioButton1->Checked==true)
	{
		rrr = false;
		delete mesh2;
		mesh2 = new Tri_Mesh;
		f_it_list.clear();
		fv.clear();
	}
	
	hkoglPanelControl1->Invalidate();
}
private: System::Void hkoglPanelControl1_MouseMove(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e)
{
	if (e->Button == System::Windows::Forms::MouseButtons::Left)
	{
		if (radioButton2->Checked == true)
		{
			point center;
			Mouse_State = Mouse::ROTATE;
			center[0] = 0.0;
			center[1] = 0.0;
			center[2] = 0.0;
			camera.mouse(e->X, e->Y, Mouse_State,
				xf * center,
				1.0, xf);
		}

		if (e->Button == System::Windows::Forms::MouseButtons::Left)
		{
			if (radioButton1->Checked == true)
			{
				boundary_rotate = 0;
#pragma region CHANGE
				GLint    viewport[4];
				GLint model[16];
				GLdouble projection[16];
				GLfloat  winX, winY, winZ;
				GLdouble posX, posY, posZ;

				glPushMatrix();

				glGetIntegerv(GL_VIEWPORT, viewport);
				glGetIntegerv(GL_MODELVIEW, model);
				glGetDoublev(GL_PROJECTION_MATRIX, projection);

				glPopMatrix();

				winX = e->X;
				winY = hkoglPanelControl1->Height - e->Y;
				glReadPixels(winX, winY, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ);
				gluUnProject(winX, winY, winZ, xf, projection, viewport, &posX, &posY, &posZ);
#pragma endregion


				if (mesh != NULL)
				{
#pragma region findPoint
					pv = mesh->vertices_begin();
					double smallest = pow(posX - mesh->point(pv).values_[0], 2) + pow(posY - mesh->point(pv).values_[1], 2) + pow(posZ - mesh->point(pv).values_[2], 2);
					for (OMT::VIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); ++v_it)
					{
						double s = pow(posX - mesh->point(v_it).values_[0], 2) + pow(posY - mesh->point(v_it).values_[1], 2) + pow(posZ - mesh->point(v_it).values_[2], 2);
						if (smallest > s)
						{
							smallest = s;
							pv = v_it;
						}
					}
#pragma endregion


					double sumf;
					double smallsumf = 99999999;

					OMT::VFIter f;
					for (OMT::VFIter vf_it = mesh->vf_iter(pv); vf_it; vf_it++)
					{
						int i = 0;
						OMT::FVIter v[3];
						double len[3];

						for (OMT::FVIter fv_it = mesh->fv_iter(vf_it); fv_it; ++fv_it)
						{
							v[i] = fv_it;

							 len[i] = TwoPointLength(posX, posY, posZ, mesh->point(fv_it).values_[0], mesh->point(fv_it).values_[1], mesh->point(fv_it).values_[2]);

							i++;
						}

						double line01, line12, line20;

						line01 = TwoPointLength(mesh->point(v[0]).values_[0], mesh->point(v[0]).values_[1], mesh->point(v[0]).values_[2], mesh->point(v[1]).values_[0], mesh->point(v[1]).values_[1], mesh->point(v[1]).values_[2]);
						line12 = TwoPointLength(mesh->point(v[1]).values_[0], mesh->point(v[1]).values_[1], mesh->point(v[1]).values_[2], mesh->point(v[2]).values_[0], mesh->point(v[2]).values_[1], mesh->point(v[2]).values_[2]);
						line20 = TwoPointLength(mesh->point(v[2]).values_[0], mesh->point(v[2]).values_[1], mesh->point(v[2]).values_[2], mesh->point(v[0]).values_[0], mesh->point(v[0]).values_[1], mesh->point(v[0]).values_[2]);

						sumf = Holen(len[0], len[1], line01);
						sumf = sumf + Holen(len[1], len[2], line12);
						sumf = sumf + +Holen(len[2], len[0], line20);
						sumf -= Holen(line01, line12, line20);

						if (sumf < smallsumf)
						{
							f = vf_it;
							smallsumf = sumf;
						}


					}

					bool flag = true;
					for (int i = 0; i < f_it_list.size(); i++)
						if (f.handle() == f_it_list[i].handle())
							flag = false;

					if (flag)
					{
						f_it_list.push_back(f);
					}
				}
			}
		}
		hkoglPanelControl1->Invalidate();
	}

	if (e->Button == System::Windows::Forms::MouseButtons::Middle)
	{
		point center;
		Mouse_State = Mouse::MOVEXY;
		center[0] = 0.0;
		center[1] = 0.0;
		center[2] = 0.0;
		camera.mouse(e->X, e->Y, Mouse_State,
			xf * center,
			1.0, xf);
		hkoglPanelControl1->Invalidate();
	}
}

#pragma region nochange



private: System::Void hkoglPanelControl1_MouseWheel(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e)
{
	if (e->Delta < 0)
	{
		point center;
		Mouse_State = Mouse::button::WHEELDOWN;
		center[0] = 0.0;
		center[1] = 0.0;
		center[2] = 0.0;
		camera.mouse(e->X, e->Y, Mouse_State,
			xf * center,
			1.0, xf);
		hkoglPanelControl1->Invalidate();
	}
	else
	{
		point center;
		Mouse_State = Mouse::WHEELUP;
		center[0] = 0.0;
		center[1] = 0.0;
		center[2] = 0.0;
		camera.mouse(e->X, e->Y, Mouse_State,
			xf * center,
			1.0, xf);
		hkoglPanelControl1->Invalidate();
	}
}
private: System::Void loadModelToolStripMenuItem_Click(System::Object^  sender, System::EventArgs^  e)
{
	openModelDialog->Filter = "Model(*.obj)|*obj";
	openModelDialog->Multiselect = false;
	openModelDialog->ShowDialog();
}
private: System::Void openModelDialog_FileOk(System::Object^  sender, System::ComponentModel::CancelEventArgs^  e)
{
	std::string filename;
	MarshalString(openModelDialog->FileName, filename);

	if (mesh != NULL)
		delete mesh;

	mesh = new Tri_Mesh;

	if (ReadFile(filename, mesh))
		std::cout << filename << std::endl;

	hkoglPanelControl1->Invalidate();
}
private: System::Void saveModelToolStripMenuItem_Click(System::Object^  sender, System::EventArgs^  e)
{
	saveModelDialog->Filter = "Model(*.obj)|*obj";
	saveModelDialog->ShowDialog();
}
private: System::Void saveModelDialog_FileOk(System::Object^  sender, System::ComponentModel::CancelEventArgs^  e)
{
	std::string filename;
	MarshalString(saveModelDialog->FileName, filename);

	if (SaveFile(filename, mesh))
		std::cout << filename << std::endl;
}
private: System::Void radioButton1_CheckedChanged(System::Object^  sender, System::EventArgs^  e) {
}
private: System::Void radioButton2_CheckedChanged(System::Object^  sender, System::EventArgs^  e) {
}

#pragma endregion

private: System::Void button1_Click(System::Object^  sender, System::EventArgs^  e)
{
	if (mesh2 != NULL&&f_it_list.size() > 0)
	{
		delete mesh2;
		mesh2 = new Tri_Mesh;
		AllPoint.clear();
		rrr = true;
		////////////////////加點////////////////////////////
		for (int i = 0; i < f_it_list.size(); i++)
		{
			for (OMT::FVIter ifv = mesh->fv_iter(f_it_list[i]); ifv; ifv++)
			{
				bool flag = true;
				for (int j = 0; j < fv.size(); j++)
				{
					if (fv[j].handle() == ifv.handle())
						flag = false;
				}

				if (flag)
					fv.push_back(ifv);
			}
		}

		for (int i = 0; i < fv.size(); i++)
		{
			mesh2->add_vertex(mesh->point(fv[i].handle()));
		}

		////////////////////加面////////////////////////////////////////////////
		for (int i = 0; i < f_it_list.size(); i++)
		{
			int ssssl = f_it_list.size();
			int sjh = ssssl;
			vector<OMT::VertexHandle> vh(3);
			int j = 0;
			for (OMT::FVIter fv_it = mesh->fv_begin(f_it_list[i]); fv_it; fv_it++)
			{
				vh[j] = fv_it;
				j++;
			}
			vector<OMT::VertexHandle> nvh(3);

			for (int k = 0; k < 3; k++)
			{
				for (OMT::VIter vi = mesh2->vertices_begin(); vi != mesh2->vertices_end(); vi++)
					if (mesh->point(vh[k]) == mesh2->point(vi))
						nvh[k] = vi.handle();
			}
			mesh2->add_face(nvh[0], nvh[1], nvh[2]);
		}

		/////////////////////讀圖////////////////////////////////
		String^s = pictureBox1->ImageLocation;
		string str;
		for (int i = 0; i < s->Length; i++)
		{
			char c = s[i];
			str += c;

		}
		cout << str;

		
		
		picturename = str;
		{
			int a=0, b=0;
			for (int i = 0; i < picturename.size(); i++)
			{
				if (picturename[i] == '\\')
				{
					a = i + 1;
					b = 0;
				}
				b++;
			}
			picturename = picturename.substr(a, b);
		}
		cout << endl << picturename << endl;
		LImage* image = loadBMP(str.c_str());
		_textureId = loadTexture(image);
		delete image;
		/////////////////////////簡易算座標/////////////////////////////////
	/*
		vector<VertexHandle> vh;
		SortBoundary(*mesh2, vh);

		DefineSidePoint(vh, *mesh2);

		vector<TexturePoint> Boundary = DefineSidePoint(vh, *mesh2);
		vector<TexturePoint> Inside = Fill_Matrix(*mesh2, Boundary);


		for (int i = 0; i < Boundary.size(); i++)
			AllPoint.push_back(Boundary[i]);
		for (int i = 0; i < Inside.size(); i++)
			AllPoint.push_back(Inside[i]);*/
			////////////////////算edge weight///////////////////////////////////

		vector<EdgeWeight> ew_list;
		for (EIter ei = mesh2->edges_begin(); ei != mesh2->edges_end(); ei++)
		{
			EdgeWeight ew;

			int face = 0;
			for (FIter fi = mesh2->faces_begin(); fi != mesh2->faces_end(); fi++)
			{
				for (FEIter fei = mesh2->fe_iter(fi); fei; fei++)
				{
					if (ei.handle() == fei.handle())
						face++;//算與edge碰到的面有幾個
				}
			}

			if (face == 2)//與edge碰到的面只有一個代表是最外圍的edge
			{
				ew.edge = ei.handle();
				ew.weight = ComputeWeight(ei.handle(), *mesh2);
				ew_list.push_back(ew);
			}
		}


		vector<VertexHandle> vh;
		SortBoundary(*mesh2, vh);

		vector<TexturePoint> Boundary = DefineSidePoint(vh, *mesh2);
		vector<TexturePoint> Inside = Advanced_Fill_Matrix(*mesh2, Boundary, ew_list);


		for (int i = 0; i < Boundary.size(); i++)
			AllPoint.push_back(Boundary[i]);
		for (int i = 0; i < Inside.size(); i++)
			AllPoint.push_back(Inside[i]);

		hkoglPanelControl1->Invalidate();
	}
}
private: System::Void button2_Click(System::Object^  sender, System::EventArgs^  e) 
{
	if (mesh2 != NULL&&rrr==true)
	{
		AllPoint.clear();
		boundary_rotate++;
		vector<EdgeWeight> ew_list;
		for (EIter ei = mesh2->edges_begin(); ei != mesh2->edges_end(); ei++)
		{
			EdgeWeight ew;

			int face = 0;
			for (FIter fi = mesh2->faces_begin(); fi != mesh2->faces_end(); fi++)
			{
				for (FEIter fei = mesh2->fe_iter(fi); fei; fei++)
				{
					if (ei.handle() == fei.handle())
						face++;//算與edge碰到的面有幾個
				}
			}

			if (face == 2)//與edge碰到的面只有一個代表是最外圍的edge
			{
				ew.edge = ei.handle();
				ew.weight = ComputeWeight(ei.handle(), *mesh2);
				ew_list.push_back(ew);
			}
		}

		vector<VertexHandle> vh;
		SortBoundary(*mesh2, vh);

		for (int i = 0; i < boundary_rotate; i++)
		{
			VertexHandle h = vh[0];
			for (int j = 0; j < vh.size() - 1; j++)
			{
				vh[j] = vh[j + 1];
			}
			vh[vh.size() - 1] = h;
		}

		vector<TexturePoint> Boundary = DefineSidePoint(vh, *mesh2);
		vector<TexturePoint> Inside = Advanced_Fill_Matrix(*mesh2, Boundary, ew_list);


		for (int i = 0; i < Boundary.size(); i++)
			AllPoint.push_back(Boundary[i]);
		for (int i = 0; i < Inside.size(); i++)
			AllPoint.push_back(Inside[i]);

		hkoglPanelControl1->Invalidate();
	}
}
private: System::Void button3_Click(System::Object^  sender, System::EventArgs^  e) 
{
	PictureLoad->Filter = "Model(*.bmp)|*bmp";
	PictureLoad->Multiselect = false;
	PictureLoad->ShowDialog();
}
private: System::Void MyForm_Load(System::Object^  sender, System::EventArgs^  e) 
{

}
private: System::Void radioButton3_CheckedChanged(System::Object^  sender, System::EventArgs^  e) 
{
	hkoglPanelControl1->Invalidate();
}
private: System::Void radioButton4_CheckedChanged(System::Object^  sender, System::EventArgs^  e) 
{
	hkoglPanelControl1->Invalidate();
}
private: System::Void PictureLoad_FileOk(System::Object^  sender, System::ComponentModel::CancelEventArgs^  e) 
{
	using namespace std;
	using namespace System;

	pictureBox1->ImageLocation = PictureLoad->FileName;

	FileData fd;
	String^ str = PictureLoad->FileName;

	string str1;
	for (int i = 0; i < str->Length; i++)////////////////////////////////////////////////
	{
		char c = str[i];
		str1 += c;
	}
	fd.filelocation=str1;
	str1 = "";
	for (int i = 0; i < str->Length; i++)///////////////////////
	{
		char c = str[i];
		str1 += str[i];
		if (c == '\\')
			str1 = "";
	}

	for (int i = 0; i < str1.size(); i++)
	{
		if (str1[i] == '.')
			str1 = str1.substr(0, i);
	}

	fd.filename = str1;
	String^ str2 = gcnew String(str1.c_str());/////////////////////////

	fd_list.push_back(fd);

	listBox1->Items->Add(str2);

}
private: System::Void listBox1_SelectedIndexChanged(System::Object^  sender, System::EventArgs^  e) 
{
	String^s = listBox1->SelectedItem->ToString();
	string str = "";
	for (int i = 0; i < s->Length; i++)
	{
		char c = s[i];
		str += c;
	}

	for(int i=0;i<fd_list.size();i++)
		if (str == fd_list[i].filename)
		{
			s= gcnew String(fd_list[i].filelocation.c_str());
			pictureBox1->ImageLocation = s;
		}

}
private: System::Void button4_Click(System::Object^  sender, System::EventArgs^  e) 
{
	if (mesh2 != NULL&&rrr==true)
	{
		TextureMesh tm;
		vector<TexturePoint> NEW;
		for (int i = 0; i < AllPoint.size(); i++)
		{
			TexturePoint tp;
			for (VIter v = mesh->vertices_begin(); v != mesh->vertices_end(); v++)
			{
				if (mesh2->point(AllPoint[i].p) == mesh->point(v.handle()))
				{
					tp.p = v.handle();
					tp.x = AllPoint[i].x;
					tp.y = AllPoint[i].y;
					break;
				}
			}
			NEW.push_back(tp);
		}
		tm.vertex = NEW;
		tm.textureID = _textureId;


		for (FIter fi = mesh2->faces_begin(); fi != mesh2->faces_end(); fi++)
		{
			Face fc;
			for (FVIter fv = mesh2->fv_iter(fi); fv; fv++)
				for (int i = 0; i < NEW.size(); i++)
					if (mesh2->point(fv.handle()) == mesh->point(NEW[i].p))
						fc.vertex.push_back(NEW[i]);

			tm.face.push_back(fc);
		}

		tm.textureName = picturename;
		setamount++;
		tm.id = setamount;
		{
			String^s = textBox1->Text;
			string name="";
			for (int i = 0; i < s->Length; i++)
			{
				char c = s[i];
				name += c;
			}
			tm.SetName = name;
		}

		if (tm.SetName == "")
		{
			tm.SetName = tm.id + 48;
			tm.SetName += "(noname)";
		}

		listBox2->Items->Add(gcnew String(tm.SetName.c_str()));

		All_texture_picture.push_back(tm);
		hkoglPanelControl1->Invalidate();
	}

}
private: System::Void SaveSet_FileOk(System::Object^  sender, System::ComponentModel::CancelEventArgs^  e) 
{
	std::string filename;
	MarshalString(SaveSet->FileName, filename);

	SaveSetFile(filename, All_texture_picture);

}
private: System::Void button5_Click(System::Object^  sender, System::EventArgs^  e) 
{
	SaveSet->Filter = "*.txt|*txt";
	SaveSet->ShowDialog();
}
private: System::Void checkBox1_CheckedChanged(System::Object^  sender, System::EventArgs^  e) 
{

}
private: System::Void listBox2_SelectedIndexChanged(System::Object^  sender, System::EventArgs^  e) 
{
	
	if (checkBox1->Checked)
	{
		Show_texture.push_back(All_texture_picture[listBox2->SelectedIndex]);
	}
	else
	{
		Show_texture.clear();
		Show_texture.push_back(All_texture_picture[listBox2->SelectedIndex]);
	}
	hkoglPanelControl1->Invalidate();
}
private: System::Void button6_Click(System::Object^  sender, System::EventArgs^  e) 
{
	LoadSet->Filter = "*.txt|*txt";
	LoadSet->Multiselect = false;
	LoadSet->ShowDialog();
}
private: System::Void LoadSet_FileOk(System::Object^  sender, System::ComponentModel::CancelEventArgs^  e) 
{
	All_texture_picture.clear();
	listBox2->Items->Clear();
	listBox1->Items->Clear();
	string filename;
	MarshalString(LoadSet->FileName, filename);
	vector<TextureMesh> load;
	LoadSetFile(filename, load);
	
	string picturelocation;
	for (int i = 0; i < load.size(); i++)
	{
		picturelocation = ".\\img\\" + load[i].textureName;
		LImage* image = loadBMP(picturelocation.c_str());
		_textureId = loadTexture(image);
		load[i].textureID = _textureId;
		delete image;


		All_texture_picture.push_back(load[i]);
		FileData fd;

		listBox2->Items->Add(gcnew String(load[i].SetName.c_str()));

		{
			string str;
			int a = 0, b = 0;
			for (int i = 0; i < picturelocation.size(); i++)
			{
				if (picturelocation[i] == '\\')
				{
					a = i + 1;
					b = 0;
				}
				b++;
			}
			str = picturelocation.substr(a, b);
			for (int i = 0; i < str.size(); i++)
				if (str[i] == '.')
				{
					str = str.substr(0, i);
					break;
				}
			fd.filename = str;
		}
		
		fd.filelocation = picturelocation;
		listBox1->Items->Add(gcnew String(fd.filename.c_str()));
		fd_list.push_back(fd);
	}


}
private: System::Void button7_Click(System::Object^  sender, System::EventArgs^  e) 
{
	Show_texture.clear();
	Show_texture = All_texture_picture;
	checkBox1->Checked = true;
	hkoglPanelControl1->Invalidate();
}
};
}


