use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum LadderKind {
    Liz500250,
    Rox400Hd,
    Gs500Rox,
}

impl LadderKind {
    pub fn display_name(self) -> &'static str {
        match self {
            Self::Liz500250 => "LIZ500_250",
            Self::Rox400Hd => "ROX400HD",
            Self::Gs500Rox => "GS500ROX",
        }
    }

    pub fn sizes(self) -> &'static [f64] {
        match self {
            Self::Liz500250 => &[
                35.0, 50.0, 75.0, 100.0, 139.0, 150.0, 160.0, 200.0, 250.0, 300.0, 340.0, 350.0,
                400.0, 450.0, 490.0, 500.0,
            ],
            Self::Rox400Hd => &[
                50.0, 60.0, 90.0, 100.0, 120.0, 150.0, 160.0, 180.0, 190.0, 200.0, 220.0, 240.0,
                260.0, 280.0, 290.0, 300.0, 320.0, 340.0, 360.0, 380.0, 400.0,
            ],
            Self::Gs500Rox => &[
                35.0, 50.0, 75.0, 100.0, 139.0, 150.0, 160.0, 200.0, 250.0, 300.0, 340.0, 350.0,
                400.0, 450.0, 490.0, 500.0,
            ],
        }
    }

    pub fn expected_peak_count(self) -> usize {
        self.sizes().len()
    }
}
