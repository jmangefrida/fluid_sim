
//use math::round;

#[derive(PartialEq)]
enum Field {
    UField,
    VField,
    SField
}

struct Fluid {
    density: f64,
    num_x: usize,
    num_y: usize,
    num_cells: i32,
    h: f64,            //size of cells
    u: Vec<f64>,      //velocity (horizontal)
    v: Vec<f64>,      //velocity (vertical)
    new_u: Vec<f64>,
    new_v: Vec<f64>,
    p: Vec<f64>,      //pressure
    s: Vec<f64>,      //scalar: 0: static object, 1: fluid
    m: Vec<f64>,      //smoke
    new_m: Vec<f64>
}

impl Fluid {
    fn new(density: f64, num_x: usize, num_y: usize, h: f64) -> Fluid {
        
        let num_cells: usize = num_x as usize * num_y as usize;
        
        Fluid { density: density,
                num_x: num_x, 
                num_y: num_y, 
                num_cells: num_cells as i32, 
                h: h, 
                u: vec![0.0; num_cells], 
                v: vec![0.0; num_cells], 
                new_u: vec![0.0; num_cells], 
                new_v: vec![0.0; num_cells], 
                p: vec![0.0; num_cells], 
                s: vec![0.0; num_cells], 
                m: vec![1.0; num_cells], 
                new_m: vec![0.0; num_cells] }

    }
    // Applys gravity to the fluid.
    // TODO:  Why skip when s is 0?
    fn integrate(&mut self, dt: f64, gravity: f64) {
        
        let n = self.num_y as usize;
        for i in 1..self.num_x as usize +1 {
            for j in 1..self.num_y as usize {
                if self.s[i*n + j] != 0.0 && self.s[i*n + j-1] != 0.0 { //if this cell, and the one below it are fluids, apply gravity
                    self.v[i*n + j] += gravity * dt;
                }
            }
        }
    }

    fn solve_incompressibility(&mut self, num_iters: i32, dt: f64) {
        let n = self.num_x;
        let cp = self.density * self.h / dt;

        for _iter in 0..num_iters {
            for i in 1..self.num_x {
                for j in 1..self.num_y {
                    if self.s[i*n + j] == 0.0 {
                        continue; //if it is a static object...
                    }

                    //let s = self.s[i*n + j];
                    let sx0 = self.s[(i-1)*n + j];
                    let sx1 = self.s[(i+1)*n + j];
                    let sy0 = self.s[i*n + j-1];
                    let sy1 = self.s[i*n + j+1];
                    let s = sx0 + sx1 + sy0 + sy1;
                    
                    if s == 0.0 {
                        continue;
                    }

                    let div = self.u[(i+1)*n + j] - self.u[(i*n +j)] + self.v[i*n + j+1] - self.v[i*n + j];
                    let p = -div / s as f64;
                    //p *= scene.over_relaxation;
                    self.p[i*n + j] += cp * p;

                    self.u[i*n + j] -= sx0 as f64 * p;
                    self.u[(i+1)*n + j] += sx1 as f64 * p;
                    self.v[i*n + j] -= sy0 as f64 * p;
                    self.v[i*n + j+1] += sy1 as f64 * p;

                }

            }
        }
    }

    fn extrapolate(&mut self) {
        let n = self.num_y;
        for i in 0..self.num_x {
            self.u[i*n] = self.u[i*n + 1];
            self.u[i*n + self.num_y - 1] = self.u[i*n + self.num_y-2];
        }
        for j in 0..self.num_y {
            self.v[j] = self.v[n + j];
            self.v[(self.num_x+1)*n + j] = self.u[(self.num_x-2)*n +j]
        }
    }

    fn sample_field(&mut self, x: f64, y: f64, field: Field) -> f64 {
        let n = self.num_y;
        let h = self.h;
        let h1 = 1.0 / h;
        let h2 = 0.5 * h;


        let x = match [x, self.num_x as f64 * h ].iter().min_by(|a, b| a.partial_cmp(b).unwrap()) {Some(val) => val.clone(), None => 0.0};
        let x: f64 = match [x, h].iter().max_by(|a,b| a.partial_cmp(b).unwrap()) {Some(val) => val.clone(), None => 0.0};

        let y = match [y, self.num_y as f64 * h ].iter().min_by(|a, b| a.partial_cmp(b).unwrap()) {Some(val) => val.clone(), None => 0.0};
        let y = match [y, h].iter().max_by(|a,b| a.partial_cmp(b).unwrap()) {Some(val) => val.clone(), None => 0.0};

        let dx: f64 = h2;
        let dy: f64 = if field == Field::SField {h2} else {0.0};
        let f = match field {Field::UField => &self.u, Field::VField => &self.v, Field::SField => &self.s};

        let x0 = match [((x - dx)*h1).floor(), (self.num_x - 1) as f64 ].iter().min_by(|a,b| a.partial_cmp(b).unwrap()) {Some(val) => val.clone(), None => 0.0};
        let tx = ((x-dx) - x0*h) *h1;
        let x1 = match [x0 + 1.0, self.num_x as f64 - 1.0].iter().min_by(|a, b| a.partial_cmp(b).unwrap()) {Some(val) => val.clone(), None => 0.0};

        let y0 = match [((y - dy)*h1).floor(), (self.num_y - 1) as f64 ].iter().min_by(|a,b| a.partial_cmp(b).unwrap()) {Some(val) => val.clone(), None => 0.0};
        let ty = ((y-dy) - y0*h) *h1;
        let y1 = match [y0 + 1.0, self.num_y as f64 - 1.0].iter().min_by(|a, b| a.partial_cmp(b).unwrap()) {Some(val) => val.clone(), None => 0.0};

        let sx = 1.0 - tx;
        let sy: f64 = 1.0 - ty;

        let val: f64 = 
        sx*sy * f[(x0 as usize * n + y0 as usize)]
        + tx*sy * f[x1 as usize*n + y0 as usize]
        + tx*ty * f[x1 as usize*n + y1 as usize]
        + sx*ty * f[x0 as usize*n + y1 as usize];

        return val;


        
        //match field {
        //    Field::UField => {let f = &mut self}
        //}


    }


    fn avg_u(&self, i: usize, j: usize) -> f64 {
        let n = self.num_y;
        let u = (self.u[(i-1)*n + j] + self.u[i*n + j] + self.u[(i-1)*n + j+1] + self.u[i*n + j+1]) * 0.25;
        return u;
    }

    fn avg_v(&self, i: usize, j:usize) -> f64 {
        let n = self.num_y;
        let v = (self.v[(i-1)*n + j] + self.v[i*n + j] + self.v[(i-1)*n + j+1] + self.v[i*n + j+1]) * 0.25;
        return v;
    }

    fn advect_vel(&mut self, dt: f64) {
        self.new_u = self.u.to_vec();
        self.new_v = self.v.to_vec();

        let n = self.num_y;
        let h = self.h;
        let h2 = 0.5 * h;

        for i in 1..self.num_x {
            for j in 1..self.num_y {
                if self.s[i*n + j] != 0.0 && self.s[(i-1)*n + j] != 0.0 && j < self.num_y -1 {
                    let mut x = i as f64 * h;
                    let mut y = j as f64 * h + h2;
                    let mut u = self.u[i*n + j];
                    let v = self.avg_v(i, j);
                    x = x - dt*u;
                    y = y - dt*v;
                    u = self.sample_field(x, y, Field::UField);
                    self.new_u[i*n + j] = u;
 
                }

                if self.s[i*n + j] != 0.0 && self.s[i*n + j-1] != 0.0 && i < self.num_x - 1 {
                    let mut x = i as f64 * h;
                    let mut y = j as f64 * h + h2;
                    let u = self.avg_u(i, j);
                    let mut v = self.v[i*n + j];
                    x = x - dt*u;
                    y = y - dt*v;
                    v = self.sample_field(x, y, Field::VField);
                    self.new_v[i*n + j] = v;
                } 
            }
        }

        self.u = self.new_u.to_vec();
        self.v = self.new_v.to_vec();

    }

    fn advect_smoke(&mut self, dt: f64) {
       self.new_m = self.m.to_vec();
       let n = self.num_y;       
       let h = self.h;
       let h2 = 0.5 * h;

       for i in 1..self.num_x {
        for j in 1..self.num_y {
            if self.s[i*n + j] != 0.0 {
                let u = (self.u[i*n +j] + self.u[(i+1)*n + j]) * 0.5;
                let v = (self.v[i*n +j] + self.v[i*n + j+1]) * 0.5;
                let x: f64 = i as f64 *h + h2 - dt*u;
                let y: f64 = i as f64 *h + h2 - dt*v;

                self.new_m[i*n + j] = self.sample_field(x, y, Field::SField);


            }
        }
       }
       self.m = self.new_m.to_vec();
    }

    fn simulate(&mut self, dt:f64, gravity: f64, num_iters: i32) {
        self.integrate(dt, gravity);
        //self.p = vec![0.0; self.num_cells];
        self.solve_incompressibility(num_iters, dt);
        self.extrapolate();
        self.advect_vel(dt);
        self.advect_smoke(dt);
    }

}