use sim::Scene;
use eframe::egui;
mod sim;

struct MyApp {
    name: String,
    age: u32,
    scene: Scene,
    frames: u64,
}

impl Default for MyApp {
    fn default() -> Self {
        Self {
            name: "Fluid Sim".to_string(),
            age: 0,
            scene: Scene {..Default::default()},
            frames: 0
        }
    }
}

impl eframe::App for MyApp {
    fn update(&mut self, ctx: &egui::Context, frame: &mut eframe::Frame) {
        
        self.scene.fluid.simulate(self.scene.dt, self.scene.gravity, self.scene.num_iters);
        self.frames += 1;

        egui::CentralPanel::default().show(ctx, |ui| {
            ui.heading("Fluid Sim");
            ui.horizontal(|ui| {
                let name_label = ui.label("text");
                ui.text_edit_singleline(&mut self.name).labelled_by(name_label.id);
            });
            ui.add(egui::Slider::new(&mut self.age, 0..=120).text("age"));
            if ui.button("Click year").clicked() {
                self.age += 1;
            }
            ui.label(format!("Hello '{}', age {}", self.name, self.frames))
        });

    }
}


fn main() {
    
    //Interface setup
    let options = eframe::NativeOptions {
        initial_window_size: Some(egui::vec2(640.0, 480.0)), ..Default::default()
    };

    eframe::run_native("Fluid Sim", options, Box::new(|_cc| Box::new(MyApp::default())));


    let mut frames: i32 = 0;
    let mut scene: Scene = Scene {..Default::default()};

    /*
    loop {
        scene.fluid.simulate(scene.dt, scene.gravity, scene.num_iters);
        frames += 1;
        println!("{}", frames);
        if frames == 60 {
            break;
        }
    }
    */
    
}

fn scene_setup(scene_nr: i32) {

}
