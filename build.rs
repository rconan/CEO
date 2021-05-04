use std::process::Command;

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    let make_ceo = Command::new("make")
        .arg("clean")
        .arg("all")
        .output()
        .expect("failed to make ceo");
    println!("status: {}", make_ceo.status);
}
