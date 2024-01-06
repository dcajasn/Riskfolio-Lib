// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2010 Gael Guennebaud <gael.guennebaud@inria.fr>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <main.h>
#include <iostream>
#include <string>

#if defined(__APPLE_CC__)
  // Prevent deprecation warnings caused by GLEW on MacOS.
  #define GL_SILENCE_DEPRECATION 1
#endif
#include <GL/glew.h>
#include <Eigen/OpenGLSupport>
#if defined(__APPLE_CC__)
  #include <GLUT/glut.h>
#else
  #include <GL/freeglut.h>
#endif

using namespace Eigen;

#define VERIFY_MATRIX(CODE,REF) { \
    glMatrixMode(GL_MODELVIEW); \
    glLoadIdentity(); \
    CODE; \
    Matrix<float,4,4,ColMajor> m; m.setZero(); \
    glGet(GL_MODELVIEW_MATRIX, m); \
    if(!(REF).cast<float>().isApprox(m)) { \
      std::cerr << "Expected:\n" << ((REF).cast<float>()) << "\n" << "got\n" << m << "\n\n"; \
    } \
    VERIFY_IS_APPROX((REF).cast<float>(), m); \
  }

#define VERIFY_UNIFORM(SUFFIX,NAME,TYPE) { \
    TYPE value; value.setRandom(); \
    TYPE data; \
    int loc = glGetUniformLocation(prg_id, #NAME); \
    VERIFY((loc!=-1) && "uniform not found"); \
    glUniform(loc,value); \
    EIGEN_CAT(glGetUniform,SUFFIX)(prg_id,loc,data.data()); \
    if(!value.isApprox(data)) { \
      std::cerr << "Expected:\n" << value << "\n" << "got\n" << data << "\n\n"; \
    } \
    VERIFY_IS_APPROX(value, data); \
  }

#define VERIFY_UNIFORMi(NAME,TYPE) { \
    TYPE value = TYPE::Random().eval().cast<float>().cast<TYPE::Scalar>(); \
    TYPE data; \
    int loc = glGetUniformLocation(prg_id, #NAME); \
    VERIFY((loc!=-1) && "uniform not found"); \
    glUniform(loc,value); \
    glGetUniformiv(prg_id,loc,(GLint*)data.data()); \
    if(!value.isApprox(data)) { \
      std::cerr << "Expected:\n" << value << "\n" << "got\n" << data << "\n\n"; \
    } \
    VERIFY_IS_APPROX(value, data); \
  }

void printProgramInfoLog(GLuint objectID)
{
    int infologLength, charsWritten;
    GLchar *infoLog;
    glGetProgramiv(objectID, GL_INFO_LOG_LENGTH, &infologLength);
    if(infologLength > 0)
    {
        infoLog = new GLchar[infologLength];
        glGetProgramInfoLog(objectID, infologLength, &charsWritten, infoLog);
        if (charsWritten > 0)
          std::cerr << "Program info : \n" << infoLog << std::endl;
        delete[] infoLog;
    }
}

void printShaderInfoLog(GLuint objectID)
{
    int infologLength, charsWritten;
    GLchar *infoLog;
    glGetShaderiv(objectID, GL_INFO_LOG_LENGTH, &infologLength);
    if(infologLength > 0)
    {
        infoLog = new GLchar[infologLength];
        glGetShaderInfoLog(objectID, infologLength, &charsWritten, infoLog);
        if (charsWritten > 0)
          std::cerr << "Shader info : \n" << infoLog << std::endl;
        delete[] infoLog;
    }
}

GLint createProgram(const char* vtx, const char* frg, bool print_errors = true)
{
  GLint prg_id = glCreateProgram();
  GLint vtx_id = glCreateShader(GL_VERTEX_SHADER);
  GLint frg_id = glCreateShader(GL_FRAGMENT_SHADER);
  GLint ok;

  glShaderSource(vtx_id, 1, &vtx, 0);
  glCompileShader(vtx_id);
  glGetShaderiv(vtx_id, GL_COMPILE_STATUS, &ok);
  if(!ok)
  {
    if (print_errors)
    {
      std::cerr << "vtx compilation failed\n";
      std::cerr << "Source:\n" << vtx << "\n";
      printShaderInfoLog(vtx_id);
    }
    glDeleteShader(vtx_id);
    return GL_ZERO;
  }

  glShaderSource(frg_id, 1, &frg, 0);
  glCompileShader(frg_id);
  glGetShaderiv(frg_id, GL_COMPILE_STATUS, &ok);
  if(!ok)
  {
    if (print_errors)
    {
      std::cerr << "frg compilation failed.\n";
      std::cerr << "Source:\n" << frg << "\n";
      printShaderInfoLog(frg_id);
    }
    glDeleteShader(vtx_id);
    glDeleteShader(frg_id);
    return GL_ZERO;
  }

  glAttachShader(prg_id, vtx_id);
  glAttachShader(prg_id, frg_id);
  glLinkProgram(prg_id);

  // Delete shaders once linked.
  glDeleteShader(vtx_id);
  glDeleteShader(frg_id);
  glGetProgramiv(prg_id, GL_LINK_STATUS, &ok);
  if(!ok)
  {
    if (print_errors)
    {
      std::cerr << "linking failed.\n";
      printProgramInfoLog(prg_id);
    }
    glDeleteProgram(prg_id);
    return GL_ZERO;
  }

  glUseProgram(prg_id);
  return prg_id;
}

GLint createProgram(const std::string& vtx, const std::string& frg, bool print_errors = true)
{
  return createProgram(vtx.c_str(), frg.c_str(), print_errors);
}

std::string getGlslVersionString(int gl_major_version, int gl_minor_version)
{
  switch (gl_major_version)
  {
    case 2:
      switch (gl_minor_version)
      {
        case 0:
          return "#version 110";
        case 1:
          return "#version 120";
      }
      break;
    case 3:
      switch (gl_minor_version)
      {
        case 0:
          return "#version 130";
        case 1:
          return "#version 140";
        case 2:
          return "#version 150";
        case 3:
          return "#version 330";
      }
      break;
    case 4:
      switch (gl_minor_version)
      {
        case 0:
          return "#version 400";
        case 1:
          return "#version 410";
        case 2:
          return "#version 420";
        case 3:
          return "#version 430";
        case 4:
          return "#version 440";
        case 5:
          return "#version 450";
        case 6:
          return "#version 460";
      }
      break;
  }
  return "";
}

void find_and_replace(
  std::string& str,
  const std::string& find,
  const std::string& replace)
{
  size_t loc = 0;
  size_t flen = find.length();
  size_t rlen = replace.length();
  while ( (loc = str.find(find, loc)) != std::string::npos) {
    str.replace(loc, flen, replace);
    loc += rlen;
  }
}

// Finds and replaces a set of substrings in a string.
std::string format(
  const std::string& str,
  const std::vector<std::string>& find,
  const std::vector<std::string>& replace)
{
  std::string out = str;
  for (std::size_t i=0; i<find.size(); ++i) {
    find_and_replace(out, find[i], replace[i]);
  }
  return out;
}

// GLUT display function that runs test.  Must be run within the display loop
// in order to properly destroy resources.
void openglsupport_test_loop()
{
  // Get context info.
  const GLubyte* gl_version_string = glGetString(GL_VERSION);
  std::cerr << "GL version: " << gl_version_string << std::endl;
  std::cerr << "GLSL version: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << std::endl;
  // Parse version from string since GL_MAJOR_VERSION is only supported in GL 3.0+.
  // Version string guaranteed to be <major>.<minor><vender extension>.
  GLint gl_major_version = gl_version_string[0] - '0';
  GLint gl_minor_version = gl_version_string[2] - '0';
  bool legacy_gl = gl_major_version < 3 || (gl_major_version == 3 && gl_minor_version < 2);

  // Fixed-function pipeline removed in OpenGL 3.2.
  if (legacy_gl)
  {
    // Draw a basic triangle.
    Vector3f v3f;
    Matrix3f rot;
    glBegin(GL_POINTS);
    {
      glVertex(v3f);
      glVertex(2*v3f+v3f);
      glVertex(rot*v3f);
    }
    glEnd();

    // 4x4 matrices
    Matrix4f mf44; mf44.setRandom();
    VERIFY_MATRIX(glLoadMatrix(mf44), mf44);
    VERIFY_MATRIX(glMultMatrix(mf44), mf44);
    Matrix4d md44; md44.setRandom();
    VERIFY_MATRIX(glLoadMatrix(md44), md44);
    VERIFY_MATRIX(glMultMatrix(md44), md44);

    // Quaternion
    Quaterniond qd(AngleAxisd(internal::random<double>(), Vector3d::Random()));
    VERIFY_MATRIX(glRotate(qd), Projective3d(qd).matrix());

    Quaternionf qf(AngleAxisf(internal::random<double>(), Vector3f::Random()));
    VERIFY_MATRIX(glRotate(qf), Projective3f(qf).matrix());

    // 3D Transform
    Transform<float,3,AffineCompact> acf3; acf3.matrix().setRandom();
    VERIFY_MATRIX(glLoadMatrix(acf3), Projective3f(acf3).matrix());
    VERIFY_MATRIX(glMultMatrix(acf3), Projective3f(acf3).matrix());

    Transform<float,3,Affine> af3(acf3);
    VERIFY_MATRIX(glLoadMatrix(af3), Projective3f(af3).matrix());
    VERIFY_MATRIX(glMultMatrix(af3), Projective3f(af3).matrix());

    Transform<float,3,Projective> pf3; pf3.matrix().setRandom();
    VERIFY_MATRIX(glLoadMatrix(pf3), Projective3f(pf3).matrix());
    VERIFY_MATRIX(glMultMatrix(pf3), Projective3f(pf3).matrix());

    Transform<double,3,AffineCompact> acd3; acd3.matrix().setRandom();
    VERIFY_MATRIX(glLoadMatrix(acd3), Projective3d(acd3).matrix());
    VERIFY_MATRIX(glMultMatrix(acd3), Projective3d(acd3).matrix());

    Transform<double,3,Affine> ad3(acd3);
    VERIFY_MATRIX(glLoadMatrix(ad3), Projective3d(ad3).matrix());
    VERIFY_MATRIX(glMultMatrix(ad3), Projective3d(ad3).matrix());

    Transform<double,3,Projective> pd3; pd3.matrix().setRandom();
    VERIFY_MATRIX(glLoadMatrix(pd3), Projective3d(pd3).matrix());
    VERIFY_MATRIX(glMultMatrix(pd3), Projective3d(pd3).matrix());

    // translations (2D and 3D)
    {
      Vector2f vf2; vf2.setRandom(); Vector3f vf23; vf23 << vf2, 0;
      VERIFY_MATRIX(glTranslate(vf2), Projective3f(Translation3f(vf23)).matrix());
      Vector2d vd2; vd2.setRandom(); Vector3d vd23; vd23 << vd2, 0;
      VERIFY_MATRIX(glTranslate(vd2), Projective3d(Translation3d(vd23)).matrix());

      Vector3f vf3; vf3.setRandom();
      VERIFY_MATRIX(glTranslate(vf3), Projective3f(Translation3f(vf3)).matrix());
      Vector3d vd3; vd3.setRandom();
      VERIFY_MATRIX(glTranslate(vd3), Projective3d(Translation3d(vd3)).matrix());

      Translation<float,3> tf3; tf3.vector().setRandom();
      VERIFY_MATRIX(glTranslate(tf3), Projective3f(tf3).matrix());

      Translation<double,3> td3;  td3.vector().setRandom();
      VERIFY_MATRIX(glTranslate(td3), Projective3d(td3).matrix());
    }

    // scaling (2D and 3D)
    {
      Vector2f vf2; vf2.setRandom(); Vector3f vf23; vf23 << vf2, 1;
      VERIFY_MATRIX(glScale(vf2), Projective3f(Scaling(vf23)).matrix());
      Vector2d vd2; vd2.setRandom(); Vector3d vd23; vd23 << vd2, 1;
      VERIFY_MATRIX(glScale(vd2), Projective3d(Scaling(vd23)).matrix());

      Vector3f vf3; vf3.setRandom();
      VERIFY_MATRIX(glScale(vf3), Projective3f(Scaling(vf3)).matrix());
      Vector3d vd3; vd3.setRandom();
      VERIFY_MATRIX(glScale(vd3), Projective3d(Scaling(vd3)).matrix());

      UniformScaling<float> usf(internal::random<float>());
      VERIFY_MATRIX(glScale(usf), Projective3f(usf).matrix());

      UniformScaling<double> usd(internal::random<double>());
      VERIFY_MATRIX(glScale(usd), Projective3d(usd).matrix());
    }
  } else {
    std::cerr << "Warning: fixed-function pipeline was not tested.\n";
  }

  // Dynamic shader substitution variables.
  // Modern shaders require a version string, and newer runtimes fail to
  // compile old GLSL versions. Thus, we dynamically set the GLSL version
  // string based on runtime. Also, pre OpenGL 3.0, the output gl_FragColor was
  // built-in. This was deprecated in OpenGL 3.0, requiring us to explicitly
  // define the output variable.
  std::vector<std::string> glsl_vars;
  glsl_vars.push_back("${GLSL_VERSION}");
  glsl_vars.push_back("${FRAG_OUTPUT_DECLARATION}");
  glsl_vars.push_back("${FRAG_OUTPUT_VARIABLE}");

  std::vector<std::string> glsl_vals;
  glsl_vals.push_back(getGlslVersionString(gl_major_version, gl_minor_version));
  if (gl_major_version >= 3) {
    glsl_vals.push_back("out vec4 fragColor;");
    glsl_vals.push_back("fragColor");
  } else {
    glsl_vals.push_back("");
    glsl_vals.push_back("gl_FragColor");
  }

  // uniform
  {
    // vertex shader.
    std::string vtx = format(
      "${GLSL_VERSION}\n"
      "void main(void) {\n"
      "  gl_Position = vec4(0,0,0,1);\n"
      "}\n",
      glsl_vars, glsl_vals);

#ifdef GL_VERSION_2_0
    if(GLEW_VERSION_2_0 && GL_VERSION_2_0)
    {
      std::string frg = format(
        "${GLSL_VERSION}\n"
        "uniform vec2 v2f;\n"
        "uniform vec3 v3f;\n"
        "uniform vec4 v4f;\n"
        "uniform ivec2 v2i;\n"
        "uniform ivec3 v3i;\n"
        "uniform ivec4 v4i;\n"
        "uniform mat2 m2f;\n"
        "uniform mat3 m3f;\n"
        "uniform mat4 m4f;\n"
        "${FRAG_OUTPUT_DECLARATION}\n"
        "void main(void) { \n"
        "  ${FRAG_OUTPUT_VARIABLE} = vec4(v2f[0]+v3f[0]+v4f[0])+vec4(v2i[0]+v3i[0]+v4i[0])+vec4(m2f[0][0]+m3f[0][0]+m4f[0][0]);\n"
        "}\n",
        glsl_vars, glsl_vals);

      GLint prg_id = createProgram(vtx, frg);
      VERIFY(prg_id > 0 && "Failed to create program.");
      VERIFY_UNIFORM(fv, v2f, Vector2f);
      VERIFY_UNIFORM(fv, v3f, Vector3f);
      VERIFY_UNIFORM(fv, v4f, Vector4f);
      VERIFY_UNIFORMi(v2i, Vector2i);
      VERIFY_UNIFORMi(v3i, Vector3i);
      VERIFY_UNIFORMi(v4i, Vector4i);
      VERIFY_UNIFORM(fv, m2f, Matrix2f);
      VERIFY_UNIFORM(fv, m3f, Matrix3f);
      VERIFY_UNIFORM(fv, m4f, Matrix4f);
      glDeleteProgram(prg_id);
    }
    else
#endif
      std::cerr << "Warning: opengl 2.0 was not tested.\n";

#ifdef GL_VERSION_2_1
    if(GLEW_VERSION_2_1 && GL_VERSION_2_1 &&
        (gl_major_version > 2 || (gl_major_version == 2 && gl_minor_version >= 1)))
    {
      std::string frg = format(
        "${GLSL_VERSION}\n"
        "uniform mat2x3 m23f;\n"
        "uniform mat3x2 m32f;\n"
        "uniform mat2x4 m24f;\n"
        "uniform mat4x2 m42f;\n"
        "uniform mat3x4 m34f;\n"
        "uniform mat4x3 m43f;\n"
        "${FRAG_OUTPUT_DECLARATION}\n"
        "void main(void) {\n"
        "  ${FRAG_OUTPUT_VARIABLE} = vec4(m23f[0][0]+m32f[0][0]+m24f[0][0]+m42f[0][0]+m34f[0][0]+m43f[0][0]);\n"
        "}\n",
        glsl_vars, glsl_vals);

      GLint prg_id = createProgram(vtx, frg);
      VERIFY(prg_id > 0 && "Failed to create program.");
      typedef Matrix<float,2,3> Matrix23f;
      typedef Matrix<float,3,2> Matrix32f;
      typedef Matrix<float,2,4> Matrix24f;
      typedef Matrix<float,4,2> Matrix42f;
      typedef Matrix<float,3,4> Matrix34f;
      typedef Matrix<float,4,3> Matrix43f;

      VERIFY_UNIFORM(fv, m23f, Matrix23f);
      VERIFY_UNIFORM(fv, m32f, Matrix32f);
      VERIFY_UNIFORM(fv, m24f, Matrix24f);
      VERIFY_UNIFORM(fv, m42f, Matrix42f);
      VERIFY_UNIFORM(fv, m34f, Matrix34f);
      VERIFY_UNIFORM(fv, m43f, Matrix43f);
      glDeleteProgram(prg_id);
    }
    else
#endif
      std::cerr << "Warning: opengl 2.1 was not tested.\n";

#ifdef GL_VERSION_3_0
    if(GLEW_VERSION_3_0 && GL_VERSION_3_0 && gl_major_version >= 3)
    {
      std::string frg = format(
        "${GLSL_VERSION}\n"
        "uniform uvec2 v2ui;\n"
        "uniform uvec3 v3ui;\n"
        "uniform uvec4 v4ui;\n"
        "${FRAG_OUTPUT_DECLARATION}\n"
        "void main(void) {\n"
        "  ${FRAG_OUTPUT_VARIABLE} = vec4(v2ui[0]+v3ui[0]+v4ui[0]);\n"
        "}\n",
        glsl_vars, glsl_vals);

      GLint prg_id = createProgram(vtx, frg);
      VERIFY(prg_id > 0 && "Failed to create program.");
      typedef Matrix<unsigned int,2,1> Vector2ui;
      typedef Matrix<unsigned int,3,1> Vector3ui;
      typedef Matrix<unsigned int,4,1> Vector4ui;

      VERIFY_UNIFORMi(v2ui, Vector2ui);
      VERIFY_UNIFORMi(v3ui, Vector3ui);
      VERIFY_UNIFORMi(v4ui, Vector4ui);
      glDeleteProgram(prg_id);
    }
    else
#endif
      std::cerr << "Warning: opengl 3.0 was not tested.\n";

    // dvecn supported if >= 4.1 or ARB_vertex_attrib_64bit
    bool has_fp64_native = (gl_major_version == 4 && gl_minor_version >= 1);
    bool has_fp64_extension = false;
#ifdef GLEW_ARB_gpu_shader_fp64
    if(GLEW_ARB_gpu_shader_fp64)
    {
      // Check that extension can actually be compiled.
      if (has_fp64_extension)
      {
        std::string frg = format(
          "${GLSL_VERSION}\n"
          "#extension GL_ARB_gpu_shader_fp64 : enable\n"
          "uniform dvec2 dv2;\n"
          "${FRAG_OUTPUT_DECLARATION}\n"
          "void main(void) {\n"
          "  ${FRAG_OUTPUT_VARIABLE} = vec4(dv2.x, dv2.y, dv2.x, dv2.y);\n"
          "}\n",
          glsl_vars, glsl_vals);
        GLint prg_id = createProgram(vtx, frg, /*print_errors=*/false);
        if (prg_id)
        {
          has_fp64_extension = true;
          glDeleteProgram(prg_id);
        }
      }
    }
#endif

    if( has_fp64_native || has_fp64_extension )
    {
      std::vector<std::string> glsl_vars_with_extension = glsl_vars;
      glsl_vars_with_extension.push_back("${GLSL_EXTENSIONS}");
      std::vector<std::string> glsl_vals_with_extension = glsl_vals;
      if (has_fp64_extension)
      {
        glsl_vals_with_extension.push_back("#extension GL_ARB_gpu_shader_fp64 : enable");
      }
      else
      {
        glsl_vals_with_extension.push_back("");
      }

      std::string frg = format(
        "${GLSL_VERSION}\n"
        "${GLSL_EXTENSIONS}\n"
        "uniform dvec2 v2d;\n"
        "uniform dvec3 v3d;\n"
        "uniform dvec4 v4d;\n"
        "${FRAG_OUTPUT_DECLARATION}\n"
        "void main(void) {\n"
        "  ${FRAG_OUTPUT_VARIABLE} = vec4(v2d[0]+v3d[0]+v4d[0]);\n"
        "}\n",
        glsl_vars_with_extension, glsl_vals_with_extension);

      GLint prg_id = createProgram(vtx,frg);
      VERIFY(prg_id > 0 && "Failed to create program.");
      VERIFY_UNIFORM(dv, v2d, Vector2d);
      VERIFY_UNIFORM(dv, v3d, Vector3d);
      VERIFY_UNIFORM(dv, v4d, Vector4d);
      glDeleteProgram(prg_id);
    }
    else
      std::cerr << "Warning: dvec (fp64) was not tested.\n";
  }

  // Exit loop - Leaving main loop is supported by freeglut, otherwise we
  // are forced to exit.
#ifdef FREEGLUT
  glutLeaveMainLoop();
  // Trigger another display loop iteration. Otherwise, it just hangs.
  glutPostRedisplay();
#else
  exit(0);
#endif
}

EIGEN_DECLARE_TEST(openglsupport)
{
  int argc = 0;
  glutInit(&argc, 0);

  GLint glut_display_mode = GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH;

#ifndef EIGEN_LEGACY_OPENGL
  // Initialize 3.2+ OpenGL context.
#if defined(__APPLE_CC__)
  glut_display_mode |= GLUT_3_2_CORE_PROFILE;
#elif defined(FREEGLUT)
  glutInitContextVersion(3, 2);
  glutInitContextFlags(GLUT_FORWARD_COMPATIBLE);
  glutInitContextProfile(GLUT_CORE_PROFILE);
#endif
#endif

  glutInitDisplayMode(glut_display_mode);
  glutInitWindowPosition(0, 0);
  glutInitWindowSize(10, 10);

  int window = glutCreateWindow("Eigen");
  if(window <= 0)
  {
    std::cerr << "Error: Unable to create GLUT Window.\n";
    exit(1);
  }

  glewExperimental = GL_TRUE;
  if(glewInit() != GLEW_OK)
  {
    std::cerr << "Warning: Failed to initialize GLEW.\n";
    exit(1);
  }

  // Run test in display, otherwise GLUT fails to clean up and leads to memory
  // access errors on exit.
  glutDisplayFunc(openglsupport_test_loop);
  glutMainLoop();
  glutDestroyWindow(window);
}
