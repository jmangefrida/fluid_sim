use sim::Scene;

mod sim;

fn main() {
    
    let mut frames: i32 = 0;
    let mut scene: Scene = Scene {..Default::default()};

    loop {
        scene.fluid.simulate(scene.dt, scene.gravity, scene.num_iters);
        frames += 1;
        println!("{}", frames);
        if frames == 60 {
            break;
        }
    }
}

fn scene_setup(scene_nr: i32) {

}
