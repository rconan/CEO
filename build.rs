use std::process::Command;

fn main() {
    let make_ceo = Command::new("make")
        .arg("all")
        .output()
        .expect("failed to make ceo");
    println!("status: {}", make_ceo.status);
}
