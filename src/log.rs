use log::LevelFilter;
use log4rs::{
    append::{console::ConsoleAppender, file::FileAppender},
    config::{Appender, Config, Root},
    encode::pattern::PatternEncoder,
};

use std::path::Path;

pub const DEFAULT_LOG_FILE: &str = "_dev/default.log";

pub fn setup_default_logging() {
    // Log only to stdout.
    let stdout = ConsoleAppender::builder().build();
    let config = Config::builder()
        .appender(Appender::builder().build("stdout", Box::new(stdout)))
        .build(Root::builder().appender("stdout").build(LevelFilter::Debug))
        .unwrap();

    let _handle = log4rs::init_config(config).unwrap();
}

pub fn setup_test_logging() {
    todo!();
}

pub fn setup_file_logging<P: AsRef<Path>>(fp: P) {
    // Log to stdout and append to file log.
    let workdir = env!("CARGO_MANIFEST_DIR");
    let workdir = Path::new(&workdir);
    let pat = "{d(%Y-%m-%d %H:%M:%S)} {l} {t} - {m}{n}";
    let pattern = Box::new(PatternEncoder::new(pat));

    let logfile = FileAppender::builder()
        .encoder(pattern.clone())
        .build(workdir.join(&fp))
        .unwrap();

    let stdout = ConsoleAppender::builder().encoder(pattern).build();

    let config = Config::builder()
        .appender(Appender::builder().build("stdout", Box::new(stdout)))
        .appender(Appender::builder().build("dev", Box::new(logfile)))
        .build(
            Root::builder()
                .appender("dev")
                .appender("stdout")
                .build(LevelFilter::Debug),
        )
        .unwrap();

    let _handle = log4rs::init_config(config).unwrap();
    log::debug!("Logging to {}", fp.as_ref().to_string_lossy());
}
