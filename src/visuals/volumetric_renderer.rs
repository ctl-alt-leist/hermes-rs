//! Volumetric point renderer with additive blending.
//!
//! Renders field density as soft, overlapping point sprites. Each point
//! is a circular blob with Gaussian falloff. Additive blending makes
//! isolated points nearly transparent while dense clusters accumulate
//! to full brightness — producing a smooth, volumetric appearance
//! without explicit alpha per object.
//!
//! Depth testing is disabled during the volumetric pass so all points
//! contribute regardless of depth order.

use kiss3d::camera::Camera;
use kiss3d::context::Context;
use kiss3d::nalgebra::{Matrix4, Point3};
use kiss3d::renderer::Renderer;
use kiss3d::resource::{
    AllocationType, BufferType, Effect, GPUVec, ShaderAttribute, ShaderUniform,
};

/// Batched volumetric point renderer with additive blending.
pub struct VolumetricRenderer {
    shader: Effect,
    pos: ShaderAttribute<Point3<f32>>,
    color: ShaderAttribute<Point3<f32>>,
    proj: ShaderUniform<Matrix4<f32>>,
    view: ShaderUniform<Matrix4<f32>>,
    point_size_uniform: ShaderUniform<f32>,
    alpha_uniform: ShaderUniform<f32>,
    falloff_uniform: ShaderUniform<f32>,
    points: GPUVec<Point3<f32>>,
    point_size: f32,
    alpha: f32,
    falloff: f32,
}

impl VolumetricRenderer {
    /// Create a new volumetric renderer.
    pub fn new(point_size: f32, alpha: f32, falloff: f32) -> Self {
        let mut shader = Effect::new_from_str(VERTEX_SRC, FRAGMENT_SRC);
        shader.use_program();

        Self {
            pos: shader.get_attrib("position").unwrap(),
            color: shader.get_attrib("color").unwrap(),
            proj: shader.get_uniform("proj").unwrap(),
            view: shader.get_uniform("view").unwrap(),
            point_size_uniform: shader.get_uniform("point_size").unwrap(),
            alpha_uniform: shader.get_uniform("blob_alpha").unwrap(),
            falloff_uniform: shader.get_uniform("blob_falloff").unwrap(),
            points: GPUVec::new(Vec::new(), BufferType::Array, AllocationType::StreamDraw),
            shader,
            point_size,
            alpha,
            falloff,
        }
    }

    /// Queue a point for drawing this frame.
    ///
    /// Points are cleared after each render. Call this every frame for
    /// each visible cell.
    pub fn draw_point(&mut self, position: Point3<f32>, color: Point3<f32>) {
        if let Some(points) = self.points.data_mut() {
            points.push(position);
            points.push(color);
        }
    }

    /// Set the screen-space point size in pixels.
    pub fn set_point_size(&mut self, size: f32) {
        self.point_size = size;
    }
}

impl Renderer for VolumetricRenderer {
    fn render(&mut self, pass: usize, camera: &mut dyn Camera) {
        if self.points.len() == 0 {
            return;
        }

        let ctxt = Context::get();

        // Enable additive blending, disable depth test.
        ctxt.enable(Context::BLEND);
        ctxt.blend_func_separate(Context::SRC_ALPHA, Context::ONE, Context::ONE, Context::ONE);
        ctxt.disable(Context::DEPTH_TEST);

        self.shader.use_program();
        self.pos.enable();
        self.color.enable();

        camera.upload(pass, &mut self.proj, &mut self.view);
        self.point_size_uniform.upload(&self.point_size);
        self.alpha_uniform.upload(&self.alpha);
        self.falloff_uniform.upload(&self.falloff);

        self.color.bind_sub_buffer(&mut self.points, 1, 1);
        self.pos.bind_sub_buffer(&mut self.points, 1, 0);

        ctxt.point_size(self.point_size);
        ctxt.draw_arrays(Context::POINTS, 0, (self.points.len() / 2) as i32);

        self.pos.disable();
        self.color.disable();

        // Restore GL state.
        ctxt.disable(Context::BLEND);
        ctxt.enable(Context::DEPTH_TEST);

        // Clear points for next frame.
        if let Some(points) = self.points.data_mut() {
            points.clear();
        }
    }
}

// ============================================================================
// Shaders
// ============================================================================

const VERTEX_SRC: &str = "#version 100
    attribute vec3 position;
    attribute vec3 color;
    varying   vec3 Color;
    uniform   mat4 proj;
    uniform   mat4 view;
    uniform   float point_size;
    void main() {
        gl_Position  = proj * view * vec4(position, 1.0);
        gl_PointSize = point_size;
        Color = color;
    }";

const FRAGMENT_SRC: &str = "#version 100
#ifdef GL_FRAGMENT_PRECISION_HIGH
   precision highp float;
#else
   precision mediump float;
#endif

    varying vec3 Color;
    uniform float blob_alpha;
    uniform float blob_falloff;
    void main() {
        vec2 coord = gl_PointCoord - vec2(0.5);
        float r2 = dot(coord, coord);
        if (r2 > 0.25) discard;

        float falloff = exp(-blob_falloff * r2);
        float alpha = falloff * blob_alpha;
        gl_FragColor = vec4(Color * alpha, alpha);
    }";
